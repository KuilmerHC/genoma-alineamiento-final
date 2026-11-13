"""Experimento C3: Smith-Waterman.

Hipótesis H3: T(m, n) = Θ(mn).
  (a) Lectura fija (m = 50) y región creciente n: pendiente log-log de T vs n ≈ 1.
  (b) Cuadrado m = n: pendiente log-log de T vs n ≈ 2.
  (c) T/(mn) aproximadamente constante en ambos casos.

Uso: python -m experimentos.exp_c3_smith_waterman
"""

from c3_smith_waterman.smith_waterman import smith_waterman
from datos.generar import generar_genoma
from experimentos.bench import (COLOR_MEDIDO, COLOR_TEORIA, grafica, guardar_csv, medir,
                                pendiente_loglog, resumen, tamanos_geometricos)


def serie(pares, repeticiones):
    filas = []
    for m, n in pares:
        a = generar_genoma(m, semilla=1)
        b = generar_genoma(n, semilla=2)
        r = resumen(medir(lambda d: smith_waterman(*d), lambda: (a, b), repeticiones))
        filas.append({"m": m, "n": n, "mn": m * n, **r,
                      "ns_por_celda": r["mediana"] / (m * n) * 1e9})
        print(f"m={m:>5} n={n:>6}  T={r['mediana']:.4f}s  T/(mn)={filas[-1]['ns_por_celda']:.1f} ns")
    return filas


def main(repeticiones=5):
    print("(a) m = 50 fijo, n creciente")
    fa = serie([(50, n) for n in tamanos_geometricos(500, 7)], repeticiones)  # n = 500 … 32000
    print("(b) m = n")
    fb = serie([(n, n) for n in tamanos_geometricos(64, 6)], repeticiones)    # n = 64 … 2048
    guardar_csv("c3_sw_m_fijo", fa)
    guardar_csv("c3_sw_cuadrado", fb)

    for filas, nombre, titulo, xlabel in (
        (fa, "c3_tiempo_m_fijo", "C3 Smith-Waterman: m = 50, n creciente", "n (largo de la región)"),
        (fb, "c3_tiempo_cuadrado", "C3 Smith-Waterman: m = n", "n (= m)"),
    ):
        ns = [f["n"] for f in filas]
        med = [f["mediana"] for f in filas]
        p = pendiente_loglog(ns, med)
        print(f"{nombre}: pendiente log-log = {p:.3f}")
        c = med[-1] / filas[-1]["mn"]
        grafica(nombre, ns, [
            {"y": med, "q1": [f["q1"] for f in filas], "q3": [f["q3"] for f in filas],
             "etiqueta": "medido (mediana, Q1–Q3)", "color": COLOR_MEDIDO},
            {"y": [c * f["mn"] for f in filas], "etiqueta": "c·m·n", "color": COLOR_TEORIA,
             "estilo": "--"},
        ], titulo, xlabel, "tiempo (s)", nota=f"pendiente = {p:.3f}")

    todas = fa + fb
    grafica("c3_razon_mn", [f["mn"] for f in fa], [
        {"y": [f["ns_por_celda"] for f in fa], "etiqueta": "m = 50 fijo",
         "color": COLOR_MEDIDO, "marcador": "o"},
    ], "C3: tiempo por celda T/(mn) con m = 50", "m·n (celdas)", "ns por celda",
        loglog=False, logx=True, y_desde_cero=True)
    return todas


if __name__ == "__main__":
    main()
