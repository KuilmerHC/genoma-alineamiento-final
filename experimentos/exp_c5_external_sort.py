"""Experimento C5: ordenamiento externo sobre DiscoSimulado. E/S contadas vs predichas.

Hipótesis H5: E/S contadas = 2·(N/B)·(1 + ⌈log_{M/B−1}(N/M)⌉). Si N es múltiplo de B la
coincidencia es exacta; si no, la diferencia es solo el último bloque incompleto.

Serie A: N entre 10⁵ y 10⁶ (5 tamaños en progresión geométrica de razón ⁴√10 ≈ 1,78),
         B = 100 registros, M = N/20 (lo máximo permitido) → 1 pasada de mezcla.
Serie B: N = 10⁶, B = 100, M decreciente → 1, 2, 3, 4 y 6 pasadas de mezcla.

Uso: python -m experimentos.exp_c5_external_sort
"""

import os
import shutil
import tempfile

from c5_external_sort.external_sort import es_predichas, ordenamiento_externo, pasadas_predichas
from comun.disco_simulado import DiscoSimulado, escribir_registros
from comun.kmers import kmers_codificados
from datos.generar import generar_genoma
from experimentos.bench import (COLOR_MEDIDO, COLOR_TEORIA, grafica, guardar_csv, medir,
                                pendiente_loglog, resumen)

B = 100


def correr(N, M, dir_trabajo, registros, repeticiones):
    """Mide el ordenamiento externo de `registros` (ya generados) con memoria M."""
    ruta = os.path.join(dir_trabajo, "entrada.bin")
    entrada = DiscoSimulado(ruta, B)
    escribir_registros(entrada, registros)  # preparar el disco no se mide ni se cuenta
    stats_ult = {}

    def preparar():
        entrada.reiniciar_contadores()
        return entrada

    def ordenar(disco):
        salida, stats = ordenamiento_externo(disco, M, dir_trabajo)
        salida.cerrar()
        os.remove(salida.ruta)
        stats_ult.update(stats)

    r = resumen(medir(ordenar, preparar, repeticiones))
    entrada.cerrar()
    s = stats_ult
    return {
        "N": N, "M": M, "B": B, "M_sobre_B": M // B, "k": s["k"],
        "corridas": s["corridas_iniciales"], "pasadas": s["pasadas"],
        "pasadas_predichas": pasadas_predichas(N, M, B),
        "lecturas": s["lecturas"], "escrituras": s["escrituras"],
        "es_contadas": s["es_contadas"], "es_predichas": s["es_predichas"],
        "error_relativo": (s["es_contadas"] - s["es_predichas"]) / s["es_predichas"],
        "pico_memoria": s["pico_memoria"], **r,
    }


def imprimir(f):
    print(f"N={f['N']:>8} M={f['M']:>6} k={f['k']:>4} corridas={f['corridas']:>5} "
          f"pasadas={f['pasadas']} E/S contadas={f['es_contadas']:>8} predichas={f['es_predichas']:>10.1f} "
          f"error={f['error_relativo']:+.4%}  pico_mem={f['pico_memoria']}  T={f['mediana']:.2f}s")


def main(repeticiones=5):
    dir_trabajo = tempfile.mkdtemp(prefix="c5_")
    try:
        registros_max = kmers_codificados(generar_genoma(10**6 + 30, semilla=42), 31)
        print("Serie A: M = N/20")
        serie_a = []
        for N in (100_000, 177_800, 316_200, 562_300, 1_000_000):  # ≈ 10⁵·(⁴√10)^i
            M = (N // 20) // B * B  # N/20 redondeado a múltiplo de B
            serie_a.append(correr(N, M, dir_trabajo, registros_max[:N], repeticiones))
            imprimir(serie_a[-1])
        print("Serie B: N = 10⁶, M variable")
        serie_b = []
        for M in (50_000, 10_000, 2_000, 1_000, 500):
            serie_b.append(correr(10**6, M, dir_trabajo, registros_max, repeticiones))
            imprimir(serie_b[-1])
    finally:
        shutil.rmtree(dir_trabajo, ignore_errors=True)

    guardar_csv("c5_serie_a", serie_a)
    guardar_csv("c5_serie_b", serie_b)
    ns = [f["N"] for f in serie_a]
    p = pendiente_loglog(ns, [f["mediana"] for f in serie_a])
    print(f"Serie A: pendiente log-log del tiempo = {p:.3f}")

    grafica("c5_es_serie_a", ns, [
        {"y": [f["es_contadas"] for f in serie_a], "etiqueta": "E/S contadas",
         "color": COLOR_MEDIDO, "marcador": "o"},
        {"y": [f["es_predichas"] for f in serie_a], "etiqueta": "E/S predichas (fórmula)",
         "color": COLOR_TEORIA, "estilo": "--", "marcador": "x"},
    ], "C5: E/S con M = N/20, B = 100", "N (registros)", "transferencias de bloque")
    grafica("c5_es_serie_b", [f["pasadas"] for f in serie_b], [
        {"y": [f["es_contadas"] for f in serie_b], "etiqueta": "E/S contadas",
         "color": COLOR_MEDIDO, "marcador": "o"},
        {"y": [f["es_predichas"] for f in serie_b], "etiqueta": "E/S predichas (fórmula)",
         "color": COLOR_TEORIA, "estilo": "--", "marcador": "x"},
    ], "C5: E/S según el número de pasadas (N = 10⁶)", "pasadas de mezcla ⌈log_(M/B−1)(N/M)⌉",
        "transferencias de bloque", loglog=False, y_desde_cero=True)
    grafica("c5_tiempo_serie_a", ns, [
        {"y": [f["mediana"] for f in serie_a], "q1": [f["q1"] for f in serie_a],
         "q3": [f["q3"] for f in serie_a], "etiqueta": "medido", "color": COLOR_MEDIDO},
    ], "C5: tiempo con M = N/20", "N (registros)", "tiempo (s)", nota=f"pendiente = {p:.3f}")
    return serie_a, serie_b


if __name__ == "__main__":
    main()
