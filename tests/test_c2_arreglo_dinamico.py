import pytest

from c2_arreglo_dinamico.arreglo_dinamico import ArregloDinamico


def test_agregar_y_leer():
    a = ArregloDinamico()
    for i in range(100):
        a.agregar(i * 2)
    assert len(a) == 100
    assert all(a[i] == i * 2 for i in range(100))


def test_indice_fuera_de_rango():
    a = ArregloDinamico()
    a.agregar(1)
    with pytest.raises(IndexError):
        a[1]


@pytest.mark.parametrize("n", [1, 10, 1000, 4097, 100_000])
def test_costo_amortizado_menor_o_igual_a_3n(n):
    a = ArregloDinamico()
    for i in range(n):
        a.agregar(i)
    assert a.costo_total() <= 3 * n


def test_copias_son_potencias_de_dos():
    # capacidad inicial 1: con 17 appends se redimensiona copiando 1+2+4+8+16 elementos
    a = ArregloDinamico()
    for i in range(17):
        a.agregar(i)
    assert a.copias == 31


@pytest.mark.parametrize("n", [1, 2, 3, 5, 1000, 4097, 100_000])
def test_saldo_de_creditos_nunca_negativo(n):
    # Prueba del método contable: con 3 créditos por append el banco nunca queda en rojo.
    a = ArregloDinamico()
    for i in range(n):
        a.agregar(i)
    assert a.saldo_minimo >= 0
    assert a.saldo == 3 * n - a.costo_total()
