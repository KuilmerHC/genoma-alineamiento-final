import itertools

import pytest

from comun.kmers import codificar_kmer, decodificar_kmer, kmers_codificados
from datos.generar import generar_genoma, generar_kmers


def test_valores_conocidos():
    assert codificar_kmer("A") == 0
    assert codificar_kmer("T") == 3
    assert codificar_kmer("AC") == 0b0001
    assert codificar_kmer("TG") == 0b1110


def test_ida_y_vuelta():
    for k in (1, 5, 31):
        s = generar_genoma(k, semilla=k)
        assert decodificar_kmer(codificar_kmer(s), k) == s


def test_orden_de_enteros_igual_a_orden_alfabetico():
    todos = ["".join(p) for p in itertools.product("ACGT", repeat=4)]  # ya en orden alfabético
    codigos = [codificar_kmer(s) for s in todos]
    assert codigos == list(range(len(todos)))


def test_ventana_deslizante_igual_a_codificar_cada_ventana():
    genoma = generar_genoma(3000, semilla=7)
    for k in (1, 8, 21, 31):
        esperado = [(codificar_kmer(s), i) for s, i in generar_kmers(genoma, k)]
        assert kmers_codificados(genoma, k) == esperado


def test_k_invalido():
    with pytest.raises(ValueError):
        kmers_codificados("ACGT", 32)
    assert kmers_codificados("ACG", 5) == []
