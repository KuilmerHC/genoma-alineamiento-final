"""Experimento C4: quicksort aleatorizado frente a merge sort (C1).

Hipótesis H4:
  (a) con claves distintas, comparaciones / (2(n+1)H_n − 4n) → 1  (≈ 2n ln n);
  (b) quicksort hace ≈ 1,39·n log₂ n comparaciones contra ≈ 0,93·n log₂ n de merge sort;
  (c) con muchas claves repetidas (k = 6: a lo sumo 4096 distintas) la partición de 3 vías
      hace MENOS comparaciones que con claves distintas;
  (d) con la entrada ya ordenada el pivote aleatorio mantiene Θ(n log n) (sin n²/2).

Uso: python -m experimentos.exp_c4_quicksort [num_tamanos]   (por defecto 10: n = 1000 … 512000)
"""

import math
import sys

from c1_mergesort import merge_sort as ms
from c4_quicksort import quicksort as qs
from comun.kmers import kmers_codificados
from datos.generar import generar_genoma
from experimentos.bench import (COLOR_EXTRA, COLOR_MEDIDO, COLOR_TEORIA, grafica,
                                guardar_csv, medir, pendiente_loglog, resumen,
                                tamanos_geometricos)


def medir_quicksort(datos_base, repeticiones):
    """Tiempo y comparaciones de quicksort; cada repetición usa otra semilla de pivotes."""
    comps = []
    semilla = [0]

    def correr(datos):
        qs.reiniciar_contador()
        qs.quicksort(datos, semilla[0])
        comps.append(qs.comparaciones)
        semilla[0] += 1

    tiempos = medir(correr, lambda: list(datos_base), repeticiones)
    return resumen(tiempos), sum(comps) / len(comps), min(comps), max(comps)


def main(num_tamanos=10, repeticiones=5):
    filas = []
    for n in tamanos_geometricos(1000, num_tamanos):
        genoma = generar_genoma(n + 30, semilla=42)
        distintos = kmers_codificados(genoma, 15)[:n]               # pares (kmer, pos): únicos
        repetidos = [x for x, _ in kmers_codificados(genoma, 6)][:n]  # ≤ 4096 kmers distintos
        ordenados = list(range(n))

        r_qs, c_qs, cmin, cmax = medir_quicksort(distintos, repeticiones)
        r_ms = resumen(medir(ms.merge_sort, lambda: list(distintos), repeticiones))
        ms.reiniciar_contador()
        ms.merge_sort(distintos)
        c_ms = ms.comparaciones
        r_rep, c_rep, _, _ = medir_quicksort(repetidos, repeticiones)
        r_ms_rep = resumen(medir(ms.merge_sort, lambda: list(repetidos), repeticiones))
        r_ord, c_ord, _, _ = medir_quicksort(ordenados, repeticiones)

        esperado = qs.comparaciones_esperadas(n)
        nlog = n * math.log2(n)
        filas.append({
            "n": n,
            "t_qs": r_qs["mediana"], "t_qs_q1": r_qs["q1"], "t_qs_q3": r_qs["q3"],
            "t_ms": r_ms["mediana"], "t_ms_q1": r_ms["q1"], "t_ms_q3": r_ms["q3"],
            "comp_qs_prom": c_qs, "comp_qs_min": cmin, "comp_qs_max": cmax,
            "comp_esperadas": esperado, "comp_qs_sobre_esperadas": c_qs / esperado,
            "comp_qs_sobre_nlogn": c_qs / nlog, "comp_ms_sobre_nlogn": c_ms / nlog,
            "t_qs_repetidos": r_rep["mediana"], "t_ms_repetidos": r_ms_rep["mediana"],
            "comp_qs_repetidos_sobre_nlogn": c_rep / nlog,
            "t_qs_ordenada": r_ord["mediana"], "comp_qs_ordenada_sobre_esperadas": c_ord / esperado,
        })
        f = filas[-1]
        print(f"n={n:>7}  qs={f['t_qs']:.3f}s ms={f['t_ms']:.3f}s  comp/E={f['comp_qs_sobre_esperadas']:.3f}"
              f"  qs/nlogn={f['comp_qs_sobre_nlogn']:.3f} ms/nlogn={f['comp_ms_sobre_nlogn']:.3f}"
              f"  repetidos qs/nlogn={f['comp_qs_repetidos_sobre_nlogn']:.3f}"
              f"  ordenada comp/E={f['comp_qs_ordenada_sobre_esperadas']:.3f}")

    ns = [f["n"] for f in filas]
    p_qs = pendiente_loglog(ns, [f["t_qs"] for f in filas])
    p_ms = pendiente_loglog(ns, [f["t_ms"] for f in filas])
    print(f"pendiente log-log: quicksort = {p_qs:.3f}, merge sort = {p_ms:.3f}")
    guardar_csv("c4_quicksort", filas)

    grafica("c4_tiempo_vs_c1", ns, [
        {"y": [f["t_qs"] for f in filas], "q1": [f["t_qs_q1"] for f in filas],
         "q3": [f["t_qs_q3"] for f in filas], "etiqueta": f"C4 quicksort (pend. {p_qs:.2f})",
         "color": COLOR_MEDIDO},
        {"y": [f["t_ms"] for f in filas], "q1": [f["t_ms_q1"] for f in filas],
         "q3": [f["t_ms_q3"] for f in filas], "etiqueta": f"C1 merge sort (pend. {p_ms:.2f})",
         "color": COLOR_TEORIA},
    ], "C4 vs C1: tiempo con claves distintas", "n (k-mers)", "tiempo (s)")
    grafica("c4_comparaciones_razon", ns, [
        {"y": [f["comp_qs_sobre_esperadas"] for f in filas],
         "etiqueta": "quicksort / (2(n+1)Hₙ − 4n)", "color": COLOR_MEDIDO, "marcador": "o"},
        {"y": [1.0] * len(ns), "etiqueta": "teoría = 1", "color": COLOR_TEORIA, "estilo": "--"},
    ], "C4: comparaciones medidas / esperadas", "n", "razón", loglog=False, logx=True,
        y_desde_cero=True)
    grafica("c4_comparaciones_vs_c1", ns, [
        {"y": [f["comp_qs_sobre_nlogn"] for f in filas], "etiqueta": "C4 quicksort, claves distintas",
         "color": COLOR_MEDIDO, "marcador": "o"},
        {"y": [f["comp_ms_sobre_nlogn"] for f in filas], "etiqueta": "C1 merge sort",
         "color": COLOR_TEORIA, "marcador": "o"},
        {"y": [f["comp_qs_repetidos_sobre_nlogn"] for f in filas],
         "etiqueta": "C4 quicksort, k-mers repetidos (k = 6)", "color": COLOR_EXTRA, "marcador": "o"},
    ], "Comparaciones / (n log₂ n)", "n", "comparaciones / (n log₂ n)", loglog=False,
        logx=True, y_desde_cero=True)
    return filas


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 10)
