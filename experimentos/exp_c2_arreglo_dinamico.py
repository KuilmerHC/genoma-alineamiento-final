"""Experimento C2: arreglo dinámico con duplicación.

Hipótesis H2: el costo acumulado de n appends (escrituras + copias) es ≤ 3n, el saldo del
método contable nunca es negativo y el tiempo total es Θ(n) (pendiente log-log ≈ 1,
tiempo por append aproximadamente constante).

Uso: python -m experimentos.exp_c2_arreglo_dinamico [num_tamanos]   (por defecto 11: hasta ~10⁶)
"""

import sys

from c2_arreglo_dinamico.arreglo_dinamico import ArregloDinamico
from experimentos.bench import (COLOR_MEDIDO, COLOR_TEORIA, grafica, guardar_csv, medir,
                                pendiente_loglog, resumen, tamanos_geometricos)


def llenar(n):
    a = ArregloDinamico()
    for i in range(n):
        a.agregar(i)
    return a


def main(num_tamanos=11, repeticiones=5):
    filas = []
    for n in tamanos_geometricos(1000, num_tamanos):
        r = resumen(medir(llenar, lambda: n, repeticiones))
        a = llenar(n)
        # Peor caso de la razón costo/n: justo después de duplicar (n = 2^j + 1).
        n_peor = 1 << n.bit_length()
        n_peor += 1
        b = llenar(n_peor)
        filas.append({"n": n, **r, "ns_por_append": r["mediana"] / n * 1e9,
                      "costo_total": a.costo_total(), "costo_sobre_n": a.costo_total() / n,
                      "saldo_minimo": a.saldo_minimo, "n_peor": n_peor,
                      "costo_sobre_n_peor": b.costo_total() / n_peor})
        print(f"n={n:>8}  T={r['mediana']:.4f}s  costo/n={a.costo_total() / n:.3f}"
              f"  saldo_min={a.saldo_minimo}  costo/n en n={n_peor}: {b.costo_total() / n_peor:.3f}")
    ns = [f["n"] for f in filas]
    p = pendiente_loglog(ns, [f["mediana"] for f in filas])
    print(f"pendiente log-log = {p:.3f}")
    guardar_csv("c2_arreglo_dinamico", filas)

    grafica("c2_tiempo_loglog", ns, [
        {"y": [f["mediana"] for f in filas], "q1": [f["q1"] for f in filas],
         "q3": [f["q3"] for f in filas], "etiqueta": "medido (mediana, Q1–Q3)",
         "color": COLOR_MEDIDO},
        {"y": [filas[-1]["mediana"] / ns[-1] * n for n in ns], "etiqueta": "c·n",
         "color": COLOR_TEORIA, "estilo": "--"},
    ], "C2 arreglo dinámico: tiempo de n appends", "n (appends)", "tiempo (s)",
        nota=f"pendiente = {p:.3f}")

    # Costo acumulado operación por operación: se ven los saltos al duplicar y que la
    # recta 3i siempre queda por encima.
    m = 300
    a = ArregloDinamico()
    acumulado, saldo = [], []
    for i in range(m):
        a.agregar(i)
        acumulado.append(a.costo_total())
        saldo.append(a.saldo)
    xs = list(range(1, m + 1))
    grafica("c2_costo_acumulado", xs, [
        {"y": acumulado, "etiqueta": "costo real acumulado", "color": COLOR_MEDIDO},
        {"y": [3 * i for i in xs], "etiqueta": "3·i (créditos cobrados)",
         "color": COLOR_TEORIA, "estilo": "--"},
    ], "C2: costo acumulado vs. créditos del método contable", "operación i",
        "costo (escrituras + copias)", loglog=False)
    grafica("c2_saldo", xs, [
        {"y": saldo, "etiqueta": "saldo", "color": COLOR_MEDIDO},
    ], "C2: saldo del banco de créditos (nunca < 0)", "operación i", "créditos",
        loglog=False)
    return filas, p


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
