"""Experimento C1: merge sort de k-mers.

Hipótesis H1: T(n) = Θ(n log n). Se espera:
  (a) pendiente log-log de T(n) ligeramente mayor que 1 (≈1.05–1.15 en este rango de n);
  (b) T(n)/(n log₂ n) aproximadamente constante al crecer n;
  (c) comparaciones contadas entre log₂(n!) (cota inferior) y n·log₂ n (cota superior).

Uso: python -m experimentos.exp_c1_mergesort [num_tamanos]   (por defecto 10: n = 1000 … 512000)
"""

import math
import sys

from c1_mergesort import merge_sort as ms
from comun.kmers import kmers_codificados
from datos.generar import generar_genoma
from experimentos.bench import (COLOR_EXTRA, COLOR_MEDIDO, COLOR_TEORIA, grafica,
                                guardar_csv, medir, pendiente_loglog, resumen,
                                tamanos_geometricos)

K = 15


def main(num_tamanos=10, repeticiones=5):
    filas = []
    for n in tamanos_geometricos(1000, num_tamanos):
        # Datos fuera del tiempo medido: n k-mers de un genoma de semilla fija.
        kmers = kmers_codificados(generar_genoma(n + K - 1, semilla=42), K)
        tiempos = medir(ms.merge_sort, lambda: list(kmers), repeticiones)
        ms.reiniciar_contador()
        ms.merge_sort(kmers)
        r = resumen(tiempos)
        nlogn = n * math.log2(n)
        filas.append({
            "n": n, **r,
            "t_sobre_nlogn_ns": r["mediana"] / nlogn * 1e9,
            "comparaciones": ms.comparaciones,
            "comp_sobre_nlogn": ms.comparaciones / nlogn,
            "log2_factorial": math.lgamma(n + 1) / math.log(2),
        })
        print(f"n={n:>8}  T={r['mediana']:.4f}s  T/(n log n)={filas[-1]['t_sobre_nlogn_ns']:.2f} ns"
              f"  comp/(n log n)={filas[-1]['comp_sobre_nlogn']:.3f}")

    ns = [f["n"] for f in filas]
    med = [f["mediana"] for f in filas]
    p = pendiente_loglog(ns, med)
    print(f"pendiente log-log = {p:.3f}")
    guardar_csv("c1_mergesort", filas)

    c = med[-1] / (ns[-1] * math.log2(ns[-1]))  # constante anclada en el n más grande
    grafica("c1_tiempo_loglog", ns, [
        {"y": med, "q1": [f["q1"] for f in filas], "q3": [f["q3"] for f in filas],
         "etiqueta": "medido (mediana, Q1–Q3)", "color": COLOR_MEDIDO},
        {"y": [c * n * math.log2(n) for n in ns], "etiqueta": "c·n log₂ n",
         "color": COLOR_TEORIA, "estilo": "--"},
    ], "C1 merge sort: tiempo", "n (k-mers)", "tiempo (s)", nota=f"pendiente = {p:.3f}")
    grafica("c1_razon_nlogn", ns, [
        {"y": [f["t_sobre_nlogn_ns"] for f in filas], "etiqueta": "T/(n log₂ n)",
         "color": COLOR_MEDIDO, "marcador": "o"},
    ], "C1 merge sort: T(n) / (n log₂ n)", "n (k-mers)", "ns por n·log₂ n",
        loglog=False, logx=True, y_desde_cero=True)
    grafica("c1_comparaciones", ns, [
        {"y": [f["comparaciones"] for f in filas], "etiqueta": "comparaciones contadas",
         "color": COLOR_MEDIDO, "marcador": "o"},
        {"y": [n * math.log2(n) for n in ns], "etiqueta": "n log₂ n (cota superior)",
         "color": COLOR_TEORIA, "estilo": "--"},
        {"y": [f["log2_factorial"] for f in filas], "etiqueta": "log₂ n! (cota inferior)",
         "color": COLOR_EXTRA, "estilo": ":"},
    ], "C1 merge sort: comparaciones", "n (k-mers)", "comparaciones")
    return filas, p


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 10)
