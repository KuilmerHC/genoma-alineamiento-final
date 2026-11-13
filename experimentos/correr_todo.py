"""Ejecuta todos los experimentos en orden, uno a la vez (para que no compitan por la CPU).

Uso: python -m experimentos.correr_todo        (≈ 20 minutos en total)
Genera experimentos/resultados/*.csv y experimentos/graficas/*.png.
"""

import importlib
import time

EXPERIMENTOS = [
    "exp_c1_mergesort",
    "exp_c2_arreglo_dinamico",
    "exp_c3_smith_waterman",
    "exp_c4_quicksort",
    "exp_c5_external_sort",
    "exp_c6_mapreduce",
    "exp_casos_y_memoria",
]


def main():
    inicio = time.perf_counter()
    for nombre in EXPERIMENTOS:
        print(f"\n===== {nombre} =====", flush=True)
        t0 = time.perf_counter()
        importlib.import_module(f"experimentos.{nombre}").main()
        print(f"({nombre}: {time.perf_counter() - t0:.0f} s)", flush=True)
    print(f"\nTodo listo en {(time.perf_counter() - inicio) / 60:.1f} min.")


if __name__ == "__main__":
    main()
