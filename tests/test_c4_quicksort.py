import random

import pytest

from c4_quicksort import quicksort as qs
from comun.kmers import kmers_codificados
from datos.generar import generar_genoma


def test_casos_borde():
    assert qs.quicksort([]) == []
    assert qs.quicksort([7]) == [7]
    assert qs.quicksort([2, 1]) == [1, 2]


@pytest.mark.parametrize("semilla", range(10))
def test_igual_que_sorted(semilla):
    rng = random.Random(semilla)
    a = [rng.randint(0, 50) for _ in range(rng.randint(0, 300))]  # con muchos repetidos
    assert qs.quicksort(list(a), semilla) == sorted(a)  # sorted solo como oráculo


def test_ordena_en_su_lugar():
    a = [3, 1, 2]
    qs.quicksort(a, 0)
    assert a == [1, 2, 3]


def test_kmers_codificados_con_duplicados():
    kmers = [x for x, _ in kmers_codificados(generar_genoma(5000, 3), 4)]  # ≤ 256 claves
    assert qs.quicksort(list(kmers), 1) == sorted(kmers)


def test_pares_kmer_posicion():
    pares = kmers_codificados(generar_genoma(3000, 4), 10)
    assert qs.quicksort(list(pares), 2) == sorted(pares)


def test_entrada_ordenada_e_invertida_no_es_cuadratica():
    n = 2000
    for a in (list(range(n)), list(range(n, 0, -1))):
        qs.reiniciar_contador()
        qs.quicksort(a, 5)
        assert a == sorted(a)
        assert qs.comparaciones < 3 * qs.comparaciones_esperadas(n)  # lejos de n²/2 = 2·10⁶


def test_todos_iguales_es_lineal():
    n = 10_000
    qs.reiniciar_contador()
    qs.quicksort([7] * n, 0)
    assert qs.comparaciones == n - 1  # una sola partición: todo cae en la zona "="


def test_promedio_de_comparaciones_cerca_del_valor_esperado():
    # Promedio sobre 30 semillas con claves distintas: debe estar a ±5 % de 2(n+1)H_n − 4n.
    n = 2000
    base = list(range(n))
    random.Random(9).shuffle(base)
    total = 0
    for s in range(30):
        qs.reiniciar_contador()
        qs.quicksort(list(base), s)
        total += qs.comparaciones
    promedio = total / 30
    esperado = qs.comparaciones_esperadas(n)
    assert abs(promedio - esperado) / esperado < 0.05


def test_formula_esperada_con_valores_chicos():
    # n = 2: siempre 1 comparación. n = 3: E = 8/3 (calculado a mano: 2(4)(11/6) − 12).
    assert qs.comparaciones_esperadas(2) == pytest.approx(1.0)
    assert qs.comparaciones_esperadas(3) == pytest.approx(8 / 3)
