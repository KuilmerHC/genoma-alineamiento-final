import math
import random

from c1_mergesort import merge_sort as ms
from comun.kmers import kmers_codificados
from datos.generar import generar_genoma, generar_kmers


def test_vacio_y_unitario():
    assert ms.merge_sort([]) == []
    assert ms.merge_sort([("A", 0)]) == [("A", 0)]


def test_ordena_igual_que_referencia():
    kmers = generar_kmers(generar_genoma(2000, 1), 8)
    assert ms.merge_sort(kmers) == sorted(kmers)  # sorted solo como oráculo


def test_ordena_kmers_codificados():
    # Caso real del índice: pares (kmer en 2 bits, posición), con muchos k-mers repetidos.
    kmers = kmers_codificados(generar_genoma(5000, 2), 6)
    assert ms.merge_sort(kmers) == sorted(kmers)


def test_ordenado_y_inverso():
    a = [(f"{i:04d}", i) for i in range(500)]
    assert ms.merge_sort(a) == a
    assert ms.merge_sort(a[::-1]) == a


def test_comparaciones_entrada_aleatoria():
    # Claves distintas en orden aleatorio (caso promedio).
    n = 1024
    a = list(range(n))
    random.Random(4).shuffle(a)
    ms.reiniciar_contador()
    ms.merge_sort(a)
    assert ms.comparaciones <= n * math.log2(n)
    # log2(n!) acota el número PROMEDIO de comparaciones de cualquier algoritmo
    # por comparaciones; en una permutación aleatoria merge sort queda por encima.
    assert ms.comparaciones >= math.log2(math.factorial(n))


def test_mejor_caso_queda_bajo_log2_factorial():
    # La cota Ω(n log n) es para el peor caso (y el promedio), no para toda entrada:
    # con la entrada ya ordenada merge sort hace (n/2)·log2 n < log2(n!) comparaciones.
    n = 1024
    ms.reiniciar_contador()
    ms.merge_sort(list(range(n)))
    assert ms.comparaciones == (n // 2) * int(math.log2(n))
    assert ms.comparaciones < math.log2(math.factorial(n))


def test_peor_caso_exacto():
    # Con la entrada construida, cada mezcla de tamaño s hace s − 1 comparaciones.
    for n in (2, 8, 64, 1024):
        entrada = ms.entrada_peor_caso(list(range(n)))
        assert sorted(entrada) == list(range(n))  # es una permutación
        ms.reiniciar_contador()
        assert ms.merge_sort(entrada) == list(range(n))
        k = int(math.log2(n))
        assert ms.comparaciones == n * k - n + 1


def test_peor_caso_ejemplo():
    assert ms.entrada_peor_caso(list(range(8))) == [0, 4, 2, 6, 1, 5, 3, 7]
