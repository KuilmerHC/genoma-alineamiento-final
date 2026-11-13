"""C5: ordenamiento externo (external merge sort) de registros (kmer, posicion).

Modelo de memoria externa: N registros en disco, memoria interna de M registros, bloques de
B registros. Una E/S = leer o escribir UN bloque (contadas por DiscoSimulado).

Fase 1, corridas iniciales: se leen M registros (M/B bloques), se ordenan en memoria con el
quicksort de C4 (en su lugar, sin memoria extra de registros) y se escriben como una
"corrida" ordenada. Quedan r₀ = ⌈N/M⌉ corridas.
    E/S de la fase 1: N/B lecturas + N/B escrituras = 2·N/B.

Fase 2, pasadas de mezcla: se mezclan k = M/B − 1 corridas a la vez. En memoria hay un
bloque de búfer por corrida (k bloques) más un bloque de salida: k + 1 = M/B bloques = M
registros. El siguiente registro de salida es el mínimo entre las cabezas de las k corridas,
y se obtiene con un min-heap propio en O(log k).
Cada pasada divide el número de corridas por k y lee y escribe cada bloque una vez:
    número de pasadas = ⌈log_k(N/M)⌉,  E/S por pasada = 2·N/B.

Total:   E/S = 2·(N/B)·(1 + ⌈log_{M/B−1}(N/M)⌉)

Las corridas son rangos contiguos de bloques, y M es múltiplo de B, así que cada corrida
empieza en un borde de bloque. Con N/B entero el conteo coincide EXACTAMENTE con la
fórmula; si no, la única diferencia es el último bloque incompleto (⌈N/B⌉ en vez de N/B).
Se alterna entre dos archivos temporales (A → B → A …): cada pasada lee de uno y escribe en
el otro.

Tiempo de CPU: Θ(N log M) en la fase 1 (esperado, por quicksort) + Θ(N log k) por pasada.
"""

import os

from c4_quicksort.quicksort import quicksort
from c5_external_sort.heap_min import MinHeap
from comun.disco_simulado import DiscoSimulado


def pasadas_predichas(N, M, B):
    """⌈log_{M/B−1}(N/M)⌉ calculado con enteros (sin errores de redondeo del logaritmo)."""
    k = M // B - 1
    corridas = -(-N // M)
    pasadas = 0
    while corridas > 1:
        corridas = -(-corridas // k)
        pasadas += 1
    return pasadas


def es_predichas(N, M, B):
    """Fórmula del modelo: 2·(N/B)·(1 + ⌈log_{M/B−1}(N/M)⌉)."""
    return 2 * (N / B) * (1 + pasadas_predichas(N, M, B))


class Memoria:
    """Contabilidad de cuántos registros hay en memoria interna (para comprobar que ≤ M)."""

    def __init__(self, M):
        self.M = M
        self.actual = 0
        self.pico = 0

    def cargar(self, n):
        self.actual += n
        if self.actual > self.pico:
            self.pico = self.actual
        if self.actual > self.M:
            raise MemoryError(f"se superó M = {self.M} registros en memoria")

    def liberar(self, n):
        self.actual -= n


def ordenamiento_externo(entrada, M, dir_trabajo):
    """Ordena los registros del disco `entrada` usando a lo sumo M registros de memoria.

    Devuelve (disco_ordenado, estadisticas). `disco_ordenado` es uno de los dos temporales
    creados en `dir_trabajo`; el llamador debe cerrarlo. Las E/S contadas son la suma de
    los contadores de los tres discos.
    """
    B = entrada.B
    if M % B != 0 or M // B < 3:
        raise ValueError("M debe ser múltiplo de B y M/B ≥ 3 (k = M/B − 1 ≥ 2)")
    N = entrada.num_registros()
    nb = entrada.num_bloques()
    mem = Memoria(M)
    disco_a = DiscoSimulado(os.path.join(dir_trabajo, "temp_a.bin"), B)
    disco_b = DiscoSimulado(os.path.join(dir_trabajo, "temp_b.bin"), B)

    # ---- Fase 1: corridas iniciales ordenadas -------------------------------------------
    bloques_por_corrida = M // B
    corridas = []  # (bloque_inicial, num_bloques): metadatos, no registros
    b = 0
    while b < nb:
        inicio = b
        buffer = []
        while b < nb and b - inicio < bloques_por_corrida:
            bloque = entrada.leer_bloque(b)
            mem.cargar(len(bloque))
            buffer.extend(bloque)
            b += 1
        quicksort(buffer, semilla=inicio)
        for j in range(0, len(buffer), B):
            disco_a.escribir_bloque(inicio + j // B, buffer[j:j + B])
        mem.liberar(len(buffer))
        corridas.append((inicio, b - inicio))

    # ---- Fase 2: pasadas de mezcla k-vías -----------------------------------------------
    k = M // B - 1
    origen, destino = disco_a, disco_b
    pasadas = 0
    while len(corridas) > 1:
        nuevas = []
        for g in range(0, len(corridas), k):
            grupo = corridas[g:g + k]
            _mezclar_grupo(origen, destino, grupo, mem)
            nuevas.append((grupo[0][0], sum(nb_r for _, nb_r in grupo)))
        corridas = nuevas
        origen, destino = destino, origen
        pasadas += 1

    stats = {
        "N": N, "M": M, "B": B, "k": k,
        "corridas_iniciales": -(-N // M) if N else 0,
        "pasadas": pasadas,
        "lecturas": entrada.lecturas + disco_a.lecturas + disco_b.lecturas,
        "escrituras": entrada.escrituras + disco_a.escrituras + disco_b.escrituras,
        "pico_memoria": mem.pico,
    }
    stats["es_contadas"] = stats["lecturas"] + stats["escrituras"]
    stats["es_predichas"] = es_predichas(N, M, B)
    # El disco que no tiene el resultado se cierra y se borra.
    destino.cerrar()
    os.remove(destino.ruta)
    return origen, stats


def _mezclar_grupo(origen, destino, grupo, mem):
    """Mezcla las corridas de `grupo` (ordenadas, contiguas en `origen`) en `destino`.

    Memoria: 1 bloque de búfer por corrida + 1 bloque de salida. Cada entrada del heap es
    (registro, r): el registro es la cabeza actual de la corrida r.
    """
    B = origen.B
    heap = MinHeap()
    buffers, pos, sig_bloque, fin_bloque = [], [], [], []
    for r, (ini, nbl) in enumerate(grupo):
        bloque = origen.leer_bloque(ini)
        mem.cargar(len(bloque))
        buffers.append(bloque)
        pos.append(1)
        sig_bloque.append(ini + 1)
        fin_bloque.append(ini + nbl)
        heap.insertar((bloque[0], r))

    salida = []
    bloque_salida = grupo[0][0]  # la corrida mezclada ocupa el mismo rango de bloques
    while len(heap) > 0:
        registro, r = heap.extraer_min()
        salida.append(registro)  # pasa del búfer de entrada al de salida: misma memoria
        if len(salida) == B:
            destino.escribir_bloque(bloque_salida, salida)
            bloque_salida += 1
            mem.liberar(len(salida))
            salida = []
        # Siguiente registro de la corrida r: del búfer, o leyendo su siguiente bloque.
        if pos[r] < len(buffers[r]):
            heap.insertar((buffers[r][pos[r]], r))
            pos[r] += 1
        elif sig_bloque[r] < fin_bloque[r]:
            bloque = origen.leer_bloque(sig_bloque[r])
            mem.cargar(len(bloque))
            sig_bloque[r] += 1
            buffers[r] = bloque
            heap.insertar((bloque[0], r))
            pos[r] = 1
    if salida:
        destino.escribir_bloque(bloque_salida, salida)
        mem.liberar(len(salida))
