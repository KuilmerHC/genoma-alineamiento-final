import random

import pytest

from c3_smith_waterman.smith_waterman import (
    puntaje_de_alineamiento,
    smith_waterman,
)
from datos.generar import generar_genoma, generar_lecturas


def global_needleman_wunsch(a, b, coincidencia=2, diferencia=-1, hueco=-2):
    """Oráculo: puntaje del mejor alineamiento GLOBAL de a con b (sin el 0 de SW)."""
    m, n = len(a), len(b)
    G = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        G[i][0] = i * hueco
    for j in range(1, n + 1):
        G[0][j] = j * hueco
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = coincidencia if a[i - 1] == b[j - 1] else diferencia
            G[i][j] = max(G[i - 1][j - 1] + s, G[i - 1][j] + hueco, G[i][j - 1] + hueco)
    return G[m][n]


def local_fuerza_bruta(a, b):
    """Definición de alineamiento local: el mejor alineamiento global entre TODAS las
    subcadenas de a y de b (incluida la vacía, que vale 0)."""
    mejor = 0
    for i in range(len(a) + 1):
        for i2 in range(i, len(a) + 1):
            for j in range(len(b) + 1):
                for j2 in range(j, len(b) + 1):
                    mejor = max(mejor, global_needleman_wunsch(a[i:i2], b[j:j2]))
    return mejor


def test_identicas():
    r = smith_waterman("ACGTACGT", "ACGTACGT")
    assert r.puntaje == 2 * 8
    assert r.alin_a == r.alin_b == "ACGTACGT"


def test_sin_bases_en_comun_da_cero():
    r = smith_waterman("AAAA", "CCCC")
    assert r.puntaje == 0
    assert r.alin_a == r.alin_b == ""


def test_cadena_vacia():
    assert smith_waterman("", "ACGT").puntaje == 0
    assert smith_waterman("ACGT", "").puntaje == 0


def test_ejemplo_clasico_con_otros_puntajes():
    # Ejemplo clásico de la literatura con +3 / −3 / −2: puntaje 13, GTT-AC contra GTTGAC.
    r = smith_waterman("TGTTACGG", "GGTTGACTA", coincidencia=3, diferencia=-3, hueco=-2)
    assert r.puntaje == 13
    assert (r.alin_a, r.alin_b) == ("GTT-AC", "GTTGAC")


def test_encuentra_la_lectura_dentro_de_la_region():
    region = "TTTTTTTTTT" + "ACGTTGCA" + "TTTTTTTTTT"
    r = smith_waterman("ACGTTGCA", region)
    assert r.puntaje == 16
    assert (r.inicio_b, r.fin_b) == (10, 18)


def test_alineamiento_con_hueco():
    # La lectura tiene una base borrada respecto de la región: SW debe abrir un hueco.
    r = smith_waterman("ACGTACGAGTACGT", "ACGTACGTAGTACGT")
    assert "-" in r.alin_a
    assert r.puntaje == 14 * 2 - 2


@pytest.mark.parametrize("semilla", range(40))
def test_igual_que_fuerza_bruta(semilla):
    rng = random.Random(semilla)
    a = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 6)))
    b = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 6)))
    assert smith_waterman(a, b).puntaje == local_fuerza_bruta(a, b)


@pytest.mark.parametrize("semilla", range(20))
def test_alineamiento_reconstruido_es_consistente(semilla):
    rng = random.Random(100 + semilla)
    a = "".join(rng.choice("ACGT") for _ in range(30))
    b = "".join(rng.choice("ACGT") for _ in range(40))
    r = smith_waterman(a, b)
    assert len(r.alin_a) == len(r.alin_b)
    assert r.alin_a.replace("-", "") == a[r.inicio_a:r.fin_a]
    assert r.alin_b.replace("-", "") == b[r.inicio_b:r.fin_b]
    assert puntaje_de_alineamiento(r.alin_a, r.alin_b) == r.puntaje
    assert r.celdas == 30 * 40


def test_ubica_lecturas_con_errores_en_su_posicion_real():
    genoma = generar_genoma(20_000, semilla=5)
    lecturas = generar_lecturas(genoma, 50, largo=50, tasa_error=0.03,
                                semilla=5, tasa_indel=0.02)
    margen = 25
    for pos, lectura in lecturas:
        ini = max(0, pos - margen)
        region = genoma[ini:pos + 50 + margen]
        r = smith_waterman(lectura, region)
        assert abs((ini + r.inicio_b) - (pos + r.inicio_a)) <= 3
        assert r.puntaje >= 50  # muy por encima de lo que da una región al azar


@pytest.mark.parametrize("semilla", range(20))
def test_dos_filas_da_el_mismo_puntaje(semilla):
    from c3_smith_waterman.smith_waterman import puntaje_dos_filas
    rng = random.Random(500 + semilla)
    a = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 40)))
    b = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 60)))
    assert puntaje_dos_filas(a, b) == smith_waterman(a, b).puntaje
