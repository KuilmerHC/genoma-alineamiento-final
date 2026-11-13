"""C6: alineamiento de lecturas en paralelo, estilo MapReduce, con multiprocessing.

Esquema:
  0. (serial)   construir el índice de semillas del genoma: C1 + C2, Θ(n log n).
  1. (serial)   partir las R lecturas en bloques y arrancar p procesos.
  2. MAP        (paralelo) cada proceso alinea las lecturas de sus bloques con el pipeline
                semillas → región candidata → Smith-Waterman (C3) y devuelve
                (id, posicion_estimada, puntaje) por lectura.
  3. REDUCE     (serial) concatenar los resultados parciales EN ORDEN de bloque y calcular
                estadísticas globales (alineadas, puntaje total).

Los procesos hijos se crean con "fork" (Linux): heredan el genoma y el índice que el padre
ya construyó, sin copiarlos ni serializarlos (copy-on-write). Solo viajan por tubería las
lecturas (ida) y los resultados (vuelta).

Análisis (n = largo del genoma, R = lecturas de largo m, región de largo m + 2·margen):
  Costo de una lectura:  c_L = Θ(m·log n)  (búsqueda binaria de sus ≈ m semillas)
                              + Θ(m·(m + 2·margen))  (Smith-Waterman)
  Trabajo:      W = Θ(n log n) + R·c_L + Θ(R)          (índice + map + reduce)
  Profundidad:  D = Θ(n log n) + c_L + Θ(R)            (con procesadores ilimitados el map
                                                         tarda lo que UNA lectura)
  Paralelismo:  W/D. Con p procesos: T_p ≈ T_serial + R·c_L/p  (Brent: T_p ≤ W/p + D).
  Comunicación: Θ(R·m) bytes de lecturas a los procesos + Θ(R) resultados de vuelta.
  Ley de Amdahl: si una fracción f del tiempo con 1 proceso es serial,
                 S(p) = T₁/T_p ≤ 1 / (f + (1 − f)/p)  y  S(∞) ≤ 1/f.
"""

import multiprocessing as mp
import time

from comun.alineador import alinear_lectura
from comun.indice import IndiceSemillas

# Estado compartido con los procesos hijos (heredado por fork, no se serializa).
_GENOMA = None
_INDICE = None


def _map_bloque(bloque):
    """MAP: alinea un bloque de lecturas [(id, lectura), ...] → [(id, pos, puntaje), ...]."""
    res = []
    for ident, lectura in bloque:
        pos, alin = alinear_lectura(lectura, _GENOMA, _INDICE)
        res.append((ident, pos, alin.puntaje if alin is not None else 0))
    return res


def reducir(parciales):
    """REDUCE: une los resultados parciales (ya en orden de bloque) y calcula estadísticas.

    Tiempo Θ(R). No hace falta ordenar: Pool.map devuelve los bloques en el orden en que se
    enviaron, y cada bloque conserva el orden de sus lecturas.
    """
    resultados = []
    alineadas = 0
    puntaje_total = 0
    for parcial in parciales:
        for ident, pos, puntaje in parcial:
            resultados.append((ident, pos, puntaje))
            if pos is not None:
                alineadas += 1
            puntaje_total += puntaje
    return resultados, {"lecturas": len(resultados), "alineadas": alineadas,
                        "puntaje_total": puntaje_total}


def partir_en_bloques(lecturas, num_bloques):
    """Parte la lista en `num_bloques` tramos contiguos de tamaño casi igual (difieren ≤ 1)."""
    R = len(lecturas)
    num_bloques = max(1, min(num_bloques, R))
    base, extra = divmod(R, num_bloques)
    bloques, ini = [], 0
    for b in range(num_bloques):
        fin = ini + base + (1 if b < extra else 0)
        bloques.append([(i, lecturas[i]) for i in range(ini, fin)])
        ini = fin
    return bloques


def alinear_mapreduce(genoma, lecturas, p, k=12, bloques_por_proceso=4, indice=None):
    """Alinea `lecturas` contra `genoma` con p procesos.

    Devuelve (resultados, estadisticas, tiempos). `tiempos` separa cada fase para medir la
    fracción serial: indice, arranque (crear procesos), map, reduce, total.
    Con p = 1 no se crean procesos: es la versión secuencial de referencia (T₁).
    Si se pasa `indice` ya construido, la fase 0 no se repite (tiempo 0).
    """
    global _GENOMA, _INDICE
    t = {}
    t0 = time.perf_counter()
    if indice is None:
        indice = IndiceSemillas(genoma, k)
    _GENOMA, _INDICE = genoma, indice
    t["indice"] = time.perf_counter() - t0

    bloques = partir_en_bloques(lecturas, p * bloques_por_proceso)
    if p == 1:
        t["arranque"] = 0.0
        t1 = time.perf_counter()
        parciales = [_map_bloque(b) for b in bloques]
        t["map"] = time.perf_counter() - t1
    else:
        t1 = time.perf_counter()
        pool = mp.get_context("fork").Pool(p)
        t["arranque"] = time.perf_counter() - t1
        t2 = time.perf_counter()
        parciales = pool.map(_map_bloque, bloques)
        t["map"] = time.perf_counter() - t2
        t3 = time.perf_counter()
        pool.close()
        pool.join()
        t["arranque"] += time.perf_counter() - t3  # cerrar procesos también es serial

    t4 = time.perf_counter()
    resultados, stats = reducir(parciales)
    t["reduce"] = time.perf_counter() - t4
    t["total"] = time.perf_counter() - t0
    return resultados, stats, t
