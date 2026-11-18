# Proyecto Final Integrador: Análisis de Algoritmos (SIST0175)

Grupo 3: alineamiento de lecturas genómicas sintéticas sobre un genoma de referencia sintético, usando k-mers como semillas.

Integrantes: Kuilmer Santiago Hernández Cruz (C3, C6), Juan José Berrio Español (C2, C5) y Santiago Serrano Cuellar (C1, C4).


## Estructura

| Carpeta | Contenido | Estado |
|---|---|---|
| `datos/generar.py` | Genoma y lecturas sintéticas con semilla fija | hecho |
| `c1_mergesort/` | C1: merge sort de k-mers con contador de comparaciones | hecho |
| `c2_arreglo_dinamico/` | C2: arreglo dinámico (contadores de escrituras y copias) | hecho |
| `c3_smith_waterman/` | C3: Smith-Waterman con reconstrucción del alineamiento | hecho |
| `c4_quicksort/` | C4: quicksort aleatorizado con partición de 3 vías | hecho |
| `comun/` | k-mers en 2 bits, índice de semillas (C1 + C2 + búsqueda binaria), alineador, `DiscoSimulado` | hecho |
| `c5_external_sort/` | C5: external merge sort con min-heap propio | hecho |
| `c6_mapreduce/` | C6: MapReduce con `multiprocessing` | hecho |
| `experimentos/` | `bench.py` (harness) + un script por componente; `resultados/` (CSV) y `graficas/` (PNG) | C1–C6 |
| `tests/` | Pruebas pytest por componente | C1–C6, comun |

## Ejecución

Requiere Python 3.10+. Todos los comandos se corren desde la raíz del repo.

    pip install -r requirements.txt

    python -m datos.generar 10000 100 42   # genera datos/genoma.txt y datos/lecturas.txt
    python -m datos.generar 10000 100 42 0.01   # igual, con 1 % de inserciones/borrados
    pytest -q                              # corre todas las pruebas
    python -m experimentos.exp_c1_mergesort         # ~1 min
    python -m experimentos.exp_c2_arreglo_dinamico  # ~15 s
    python -m experimentos.exp_c3_smith_waterman    # ~30 s
    python -m experimentos.exp_c4_quicksort         # ~2 min
    python -m experimentos.exp_c5_external_sort     # ~5 min
    python -m experimentos.exp_c6_mapreduce         # ~5 min (usa todos los núcleos)
    python -m experimentos.exp_casos_y_memoria      # ~3 min: mejor/peor caso de C1 y memoria
    python -m experimentos.correr_todo              # todos los anteriores, en orden (~20 min)


## Declaración de uso de IA

Se usó un asistente de IA (Claude) para apoyar el diseño, la organización del repositorio y la depuración. Todo el código fue revisado por el equipo.
