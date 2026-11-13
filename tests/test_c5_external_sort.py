import os
import random

import pytest

from c5_external_sort.external_sort import es_predichas, ordenamiento_externo, pasadas_predichas
from c5_external_sort.heap_min import MinHeap
from comun.disco_simulado import DiscoSimulado, escribir_registros, leer_todo
from comun.kmers import kmers_codificados
from datos.generar import generar_genoma


# ---------------- DiscoSimulado ----------------
def test_disco_lee_lo_que_escribe_y_cuenta(tmp_path):
    with DiscoSimulado(tmp_path / "d.bin", B=4) as d:
        regs = [(2**63 + i, i) for i in range(10)]  # kmers grandes: prueba uint64
        escribir_registros(d, regs)
        assert d.escrituras == 3  # ⌈10/4⌉ bloques
        assert d.num_registros() == 10 and d.num_bloques() == 3
        assert d.leer_bloque(2) == regs[8:]  # último bloque incompleto
        assert leer_todo(d) == regs
        assert d.lecturas == 4


def test_disco_rechaza_bloque_demasiado_grande(tmp_path):
    with DiscoSimulado(tmp_path / "d.bin", B=2) as d:
        with pytest.raises(ValueError):
            d.escribir_bloque(0, [(1, 1)] * 3)


# ---------------- MinHeap ----------------
def test_heap_extrae_en_orden():
    rng = random.Random(1)
    xs = [rng.randint(0, 100) for _ in range(500)]
    h = MinHeap()
    for x in xs:
        h.insertar(x)
    salida = [h.extraer_min() for _ in range(len(xs))]
    assert salida == sorted(xs)
    assert len(h) == 0


# ---------------- Ordenamiento externo ----------------
def _ordenar(tmp_path, registros, B, M):
    entrada = DiscoSimulado(os.path.join(tmp_path, "entrada.bin"), B)
    escribir_registros(entrada, registros)
    entrada.reiniciar_contadores()  # preparar los datos no cuenta
    salida, stats = ordenamiento_externo(entrada, M, str(tmp_path))
    contadas = stats["es_contadas"]
    resultado = leer_todo(salida)
    salida.cerrar()
    entrada.cerrar()
    return resultado, stats, contadas


@pytest.mark.parametrize("N,B,M", [
    (1000, 8, 32),     # k = 3, 32 corridas → 4 pasadas
    (4096, 16, 64),    # k = 3, 64 corridas → 4 pasadas (N/B entero)
    (5000, 10, 1000),  # k = 99, 5 corridas → 1 pasada
    (300, 10, 1000),   # cabe en memoria: 0 pasadas
    (12345, 25, 100),  # N/B no entero
])
def test_ordena_y_respeta_memoria(tmp_path, N, B, M):
    regs = kmers_codificados(generar_genoma(N + 11, semilla=N), 12)
    random.Random(N).shuffle(regs)
    resultado, stats, _ = _ordenar(tmp_path, regs, B, M)
    assert resultado == sorted(regs)  # sorted solo como oráculo
    assert stats["pico_memoria"] <= M
    assert stats["pasadas"] == pasadas_predichas(N, M, B)


@pytest.mark.parametrize("N,B,M", [(4096, 16, 64), (8192, 32, 256), (6400, 64, 640)])
def test_es_contadas_igual_a_formula_con_N_multiplo_de_B(tmp_path, N, B, M):
    regs = kmers_codificados(generar_genoma(N + 11, semilla=1), 12)
    random.Random(2).shuffle(regs)
    _, stats, contadas = _ordenar(tmp_path, regs, B, M)
    assert contadas == es_predichas(N, M, B)


def test_es_con_N_no_multiplo_de_B_solo_difiere_por_el_ultimo_bloque(tmp_path):
    N, B, M = 12345, 25, 100
    regs = kmers_codificados(generar_genoma(N + 11, semilla=3), 12)
    _, stats, contadas = _ordenar(tmp_path, regs, B, M)
    nb = -(-N // B)  # ⌈N/B⌉
    assert contadas == 2 * nb * (1 + stats["pasadas"])
    assert abs(contadas - es_predichas(N, M, B)) / contadas < 0.01


def test_duplicados_y_vacio(tmp_path):
    regs = [(7, i % 3) for i in range(200)]
    resultado, _, _ = _ordenar(tmp_path, regs, 4, 16)
    assert resultado == sorted(regs)
    resultado, stats, contadas = _ordenar(tmp_path, [], 4, 16)
    assert resultado == [] and contadas == 0


def test_pasadas_predichas():
    assert pasadas_predichas(N=1000, M=1000, B=10) == 0
    assert pasadas_predichas(N=10**6, M=50_000, B=64) == 1        # 20 corridas, k = 780
    assert pasadas_predichas(N=10**6, M=10_000, B=1000) == 3      # 100 corridas, k = 9: 100 → 12 → 2 → 1
