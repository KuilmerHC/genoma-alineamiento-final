"""Pipeline de alineamiento de una lectura: semillas → región candidata → Smith-Waterman.

Para una lectura de largo m:
  1. Cada k-mer de la lectura (desplazamiento o) se busca en el índice. Si aparece en la
     posición p del genoma, "vota" por que la lectura empieza en p − o.
  2. La posición con más votos es la candidata. Se recorta la región
     genoma[candidata − margen : candidata + m + margen]; el margen absorbe los indels.
  3. Smith-Waterman (C3) alinea la lectura contra esa región: Θ(m·(m + 2·margen)),
     en vez de Θ(m·n) contra el genoma completo. Esa es la razón de usar semillas.
"""

from c3_smith_waterman.smith_waterman import smith_waterman
from comun.kmers import kmers_codificados


def alinear_lectura(lectura, genoma, indice, margen=10):
    """Devuelve (posicion_estimada, Alineamiento) o (None, None) si ninguna semilla aparece."""
    votos = {}
    for kmer, o in kmers_codificados(lectura, indice.k):
        for p in indice.buscar(kmer):
            inicio = p - o
            votos[inicio] = votos.get(inicio, 0) + 1
    if not votos:
        return None, None

    candidata, mas_votos = None, -1
    for inicio, v in votos.items():  # máximo lineal, sin max()/sorted()
        if v > mas_votos:
            candidata, mas_votos = inicio, v

    ini = max(0, candidata - margen)
    region = genoma[ini:candidata + len(lectura) + margen]
    alin = smith_waterman(lectura, region)
    # Posición en el genoma donde empezaría la lectura completa según el alineamiento.
    return ini + alin.inicio_b - alin.inicio_a, alin
