"""Experimentos complementarios para la entrega final.

Parte A, C1 hasta n = 2²⁰ (> 10⁶) con las tres entradas que provocan cada caso:
  H1-casos: comparaciones del mejor caso (entrada ordenada) = (n/2)·log₂ n,
            del peor caso (entrada "intercalada") = n·log₂ n − n + 1, ambas EXACTAS;
            la entrada aleatoria queda entre las dos y por encima de log₂ n!;
            el TIEMPO es Θ(n log n) en los tres casos (pendientes ≈ 1,1).
Parte B, memoria extra con tracemalloc (pico de bytes reservados durante el algoritmo):
  H-mem1: C1 merge sort usa Θ(n) memoria extra; C4 quicksort (en su lugar) usa O(log n),
          es decir, prácticamente nada que crezca con n.
  H-mem3: C3 con tabla completa usa Θ(mn); con dos filas, Θ(n).

Uso: python -m experimentos.exp_casos_y_memoria
"""

import math
import random
import tracemalloc

from c1_mergesort import merge_sort as ms
from c3_smith_waterman.smith_waterman import puntaje_dos_filas, smith_waterman
from c4_quicksort.quicksort import quicksort
from datos.generar import generar_genoma
from experimentos.bench import (COLOR_EXTRA, COLOR_MEDIDO, COLOR_TEORIA, grafica, guardar_csv,
                                medir, pendiente_loglog, resumen, tamanos_geometricos)


def comparaciones(entrada):
    ms.reiniciar_contador()
    ms.merge_sort(entrada)
    return ms.comparaciones


def parte_a(repeticiones):
    filas = []
    for n in tamanos_geometricos(1024, 11):  # 2¹⁰ … 2²⁰ = 1 048 576
        ordenada = list(range(n))
        peor = ms.entrada_peor_caso(ordenada)
        aleatoria = list(range(n))
        random.Random(n).shuffle(aleatoria)
        k = int(math.log2(n))
        fila = {"n": n,
                "comp_mejor": comparaciones(ordenada), "formula_mejor": (n // 2) * k,
                "comp_peor": comparaciones(peor), "formula_peor": n * k - n + 1,
                "comp_aleatoria": comparaciones(aleatoria),
                "log2_factorial": math.lgamma(n + 1) / math.log(2)}
        for nombre, datos in (("mejor", ordenada), ("peor", peor), ("aleatoria", aleatoria)):
            r = resumen(medir(ms.merge_sort, lambda d=datos: list(d), repeticiones))
            fila[f"t_{nombre}"] = r["mediana"]
            fila[f"t_{nombre}_q1"] = r["q1"]
            fila[f"t_{nombre}_q3"] = r["q3"]
        filas.append(fila)
        print(f"n={n:>8}  mejor {fila['comp_mejor']} (= {fila['formula_mejor']})  "
              f"peor {fila['comp_peor']} (= {fila['formula_peor']})  aleatoria {fila['comp_aleatoria']}  "
              f"T mejor/peor/aleat = {fila['t_mejor']:.2f}/{fila['t_peor']:.2f}/{fila['t_aleatoria']:.2f}s")
    ns = [f["n"] for f in filas]
    pend = {c: pendiente_loglog(ns, [f[f"t_{c}"] for f in filas]) for c in ("mejor", "peor", "aleatoria")}
    print("pendientes log-log del tiempo:", {c: round(p, 3) for c, p in pend.items()})
    guardar_csv("c1_casos", filas)
    grafica("c1_casos_comparaciones", ns, [
        {"y": [f["comp_peor"] / (n * math.log2(n)) for f, n in zip(filas, ns)],
         "etiqueta": "peor caso (intercalada)", "color": COLOR_TEORIA, "marcador": "o"},
        {"y": [f["comp_aleatoria"] / (n * math.log2(n)) for f, n in zip(filas, ns)],
         "etiqueta": "aleatoria", "color": COLOR_MEDIDO, "marcador": "o"},
        {"y": [f["comp_mejor"] / (n * math.log2(n)) for f, n in zip(filas, ns)],
         "etiqueta": "mejor caso (ordenada)", "color": COLOR_EXTRA, "marcador": "o"},
    ], "C1: comparaciones / (n log₂ n) según la entrada", "n", "comparaciones / (n log₂ n)",
        loglog=False, logx=True, y_desde_cero=True)
    grafica("c1_casos_tiempo", ns, [
        {"y": [f["t_peor"] for f in filas], "etiqueta": f"peor (pend. {pend['peor']:.2f})",
         "color": COLOR_TEORIA, "marcador": "o"},
        {"y": [f["t_aleatoria"] for f in filas], "etiqueta": f"aleatoria (pend. {pend['aleatoria']:.2f})",
         "color": COLOR_MEDIDO, "marcador": "o"},
        {"y": [f["t_mejor"] for f in filas], "etiqueta": f"mejor (pend. {pend['mejor']:.2f})",
         "color": COLOR_EXTRA, "marcador": "o"},
    ], "C1: tiempo en los tres casos (n hasta 2²⁰)", "n", "tiempo (s)")
    return filas, pend


def pico_memoria(funcion, datos):
    """Bytes del pico de memoria reservada DURANTE funcion(datos) (los datos ya existían)."""
    tracemalloc.start()
    tracemalloc.reset_peak()
    base = tracemalloc.get_traced_memory()[0]
    funcion(datos)
    pico = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    return pico - base


def parte_b():
    filas_orden = []
    for n in tamanos_geometricos(1000, 9):  # 10³ … 2,56·10⁵
        datos = list(range(n))
        random.Random(1).shuffle(datos)
        m_ms = pico_memoria(ms.merge_sort, list(datos))
        m_qs = pico_memoria(lambda d: quicksort(d, 0), list(datos))
        filas_orden.append({"n": n, "bytes_merge_sort": m_ms, "bytes_quicksort": m_qs,
                            "bytes_por_elem_ms": m_ms / n, "bytes_por_elem_qs": m_qs / n})
        print(f"n={n:>7}  merge sort {m_ms / n:6.1f} B/elem   quicksort {m_qs:>6} B total")
    guardar_csv("memoria_c1_c4", filas_orden)
    ns = [f["n"] for f in filas_orden]
    p_ms = pendiente_loglog(ns, [f["bytes_merge_sort"] for f in filas_orden])
    print(f"pendiente memoria merge sort = {p_ms:.3f}")
    grafica("memoria_c1_c4", ns, [
        {"y": [f["bytes_merge_sort"] for f in filas_orden], "etiqueta": f"C1 merge sort (pend. {p_ms:.2f})",
         "color": COLOR_MEDIDO, "marcador": "o"},
        {"y": [max(f["bytes_quicksort"], 1) for f in filas_orden], "etiqueta": "C4 quicksort (en su lugar)",
         "color": COLOR_TEORIA, "marcador": "o"},
    ], "Memoria extra al ordenar (pico, tracemalloc)", "n", "bytes")

    filas_sw = []
    lectura = generar_genoma(50, semilla=1)
    for n in tamanos_geometricos(500, 7):  # región de 500 a 32 000 bases
        region = generar_genoma(n, semilla=2)
        m_tabla = pico_memoria(lambda r: smith_waterman(lectura, r), region)
        m_filas = pico_memoria(lambda r: puntaje_dos_filas(lectura, r), region)
        filas_sw.append({"m": 50, "n": n, "bytes_tabla": m_tabla, "bytes_dos_filas": m_filas,
                         "bytes_por_celda_tabla": m_tabla / (50 * n),
                         "bytes_por_columna_dos_filas": m_filas / n})
        print(f"m=50 n={n:>6}  tabla {m_tabla / 1e6:7.2f} MB   dos filas {m_filas / 1e3:7.1f} KB")
    guardar_csv("memoria_c3", filas_sw)
    ns = [f["n"] for f in filas_sw]
    grafica("memoria_c3", ns, [
        {"y": [f["bytes_tabla"] for f in filas_sw], "etiqueta": "tabla completa Θ(mn)",
         "color": COLOR_MEDIDO, "marcador": "o"},
        {"y": [f["bytes_dos_filas"] for f in filas_sw], "etiqueta": "dos filas Θ(n)",
         "color": COLOR_TEORIA, "marcador": "o"},
    ], "C3: memoria de Smith-Waterman (m = 50)", "n (largo de la región)", "bytes")
    return filas_orden, filas_sw


def main(repeticiones=5):
    parte_a(repeticiones)
    parte_b()


if __name__ == "__main__":
    main()
