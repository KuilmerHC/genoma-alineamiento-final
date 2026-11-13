from comun.alineador import alinear_lectura
from comun.indice import IndiceSemillas
from comun.kmers import codificar_kmer
from datos.generar import generar_genoma, generar_lecturas


def test_buscar_devuelve_todas_las_posiciones():
    genoma = generar_genoma(5000, semilla=3)
    k = 6
    indice = IndiceSemillas(genoma, k)
    # Oráculo: búsqueda ingenua de cada k-mer en el genoma.
    for s in {genoma[i:i + k] for i in range(0, 4000, 97)}:
        esperado = [i for i in range(len(genoma) - k + 1) if genoma[i:i + k] == s]
        assert list(indice.buscar(codificar_kmer(s))) == esperado


def test_claves_ordenadas_y_sin_repetir():
    indice = IndiceSemillas(generar_genoma(3000, semilla=4), 5)
    c = indice.claves
    assert all(c[i] < c[i + 1] for i in range(len(c) - 1))


def test_kmer_ausente():
    indice = IndiceSemillas("AAAAAAAAAA", 3)
    assert list(indice.buscar(codificar_kmer("CGT"))) == []


def test_pipeline_ubica_lecturas_con_errores():
    genoma = generar_genoma(30_000, semilla=8)
    lecturas = generar_lecturas(genoma, 100, largo=50, tasa_error=0.02,
                                semilla=8, tasa_indel=0.01)
    indice = IndiceSemillas(genoma, k=12)
    aciertos = 0
    for pos, lectura in lecturas:
        estimada, _ = alinear_lectura(lectura, genoma, indice)
        if estimada is not None and abs(estimada - pos) <= 3:
            aciertos += 1
    assert aciertos >= 95
