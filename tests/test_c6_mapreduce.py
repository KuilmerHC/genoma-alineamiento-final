import pytest

from c6_mapreduce.mapreduce import alinear_mapreduce, partir_en_bloques, reducir
from datos.generar import generar_genoma, generar_lecturas


@pytest.fixture(scope="module")
def datos():
    genoma = generar_genoma(20_000, semilla=11)
    lecturas_con_pos = generar_lecturas(genoma, 120, tasa_error=0.02, semilla=11, tasa_indel=0.01)
    return genoma, [l for _, l in lecturas_con_pos], [p for p, _ in lecturas_con_pos]


def test_partir_en_bloques_cubre_todo_en_orden():
    xs = list("abcdefghij")
    bloques = partir_en_bloques(xs, 3)
    assert [len(b) for b in bloques] == [4, 3, 3]
    assert [x for b in bloques for _, x in b] == xs
    assert len(partir_en_bloques(xs, 50)) == 10  # nunca más bloques que lecturas


def test_reduce_conserva_orden_y_cuenta():
    res, stats = reducir([[(0, 5, 10), (1, None, 0)], [(2, 7, 20)]])
    assert [r[0] for r in res] == [0, 1, 2]
    assert stats == {"lecturas": 3, "alineadas": 2, "puntaje_total": 30}


@pytest.mark.parametrize("p", [2, 4])
def test_paralelo_igual_que_secuencial(datos, p):
    genoma, lecturas, _ = datos
    sec, stats_sec, _ = alinear_mapreduce(genoma, lecturas, p=1)
    par, stats_par, t = alinear_mapreduce(genoma, lecturas, p=p)
    assert par == sec
    assert stats_par == stats_sec
    assert set(t) == {"indice", "arranque", "map", "reduce", "total"}


def test_ubica_las_lecturas(datos):
    genoma, lecturas, posiciones = datos
    res, stats, _ = alinear_mapreduce(genoma, lecturas, p=2)
    aciertos = sum(1 for (i, pos, _), real in zip(res, posiciones)
                   if pos is not None and abs(pos - real) <= 3)
    assert aciertos >= 0.95 * len(lecturas)
    assert stats["lecturas"] == len(lecturas)
