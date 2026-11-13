"""Generador de datos sintéticos (Grupo 3: alineamiento de lecturas genómicas).

Uso (desde la raíz del repo):
    python -m datos.generar <largo_genoma> <num_lecturas> [semilla] [tasa_indel]
Genera datos/genoma.txt y datos/lecturas.txt (junto a este script).
"""
import os
import random
import sys

DIR_DATOS = os.path.dirname(os.path.abspath(__file__))

ALFABETO = "ACGT"


def generar_genoma(n, semilla=42):
    rng = random.Random(semilla)
    return "".join(rng.choice(ALFABETO) for _ in range(n))


def generar_lecturas(genoma, num_lecturas, largo=50, tasa_error=0.02, semilla=42,
                     tasa_indel=0.0):
    """Lecturas de largo `largo` tomadas del genoma, con errores aleatorios.
    Devuelve una lista de (posicion_real, lectura).

    Errores simulados (como en un secuenciador real):
      - sustitución: con probabilidad `tasa_error` una base se cambia por una al azar;
      - inserción / borrado (indel): con probabilidad `tasa_indel` por paso se inserta una
        base al azar en la lectura o se salta una base del genoma (50 % cada uno).
    Con tasa_indel=0 se generan exactamente los mismos datos que la versión sin indels.
    """
    rng = random.Random(semilla + 1)
    # Con borrados la lectura consume más genoma: se deja un margen de `largo` bases.
    margen = largo if tasa_indel > 0 else 0
    lecturas = []
    for _ in range(num_lecturas):
        pos = rng.randrange(len(genoma) - largo - margen + 1)
        lectura = []
        j = pos  # posición actual en el genoma
        while len(lectura) < largo and j < len(genoma):
            if tasa_indel > 0 and rng.random() < tasa_indel:
                if rng.random() < 0.5:
                    lectura.append(rng.choice(ALFABETO))  # inserción
                else:
                    j += 1                                # borrado
                continue
            base = genoma[j]
            if rng.random() < tasa_error:
                base = rng.choice(ALFABETO)               # sustitución
            lectura.append(base)
            j += 1
        lecturas.append((pos, "".join(lectura)))
    return lecturas


def generar_kmers(texto, k):
    """Lista de (kmer, posicion) para todas las ventanas de tamaño k."""
    return [(texto[i:i + k], i) for i in range(len(texto) - k + 1)]


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000
    m = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    semilla = int(sys.argv[3]) if len(sys.argv) > 3 else 42
    tasa_indel = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    genoma = generar_genoma(n, semilla)
    lecturas = generar_lecturas(genoma, m, semilla=semilla, tasa_indel=tasa_indel)
    with open(os.path.join(DIR_DATOS, "genoma.txt"), "w") as f:
        f.write(genoma)
    with open(os.path.join(DIR_DATOS, "lecturas.txt"), "w") as f:
        for pos, lec in lecturas:
            f.write(f"{pos}\t{lec}\n")
    print(f"genoma de {n} bases y {m} lecturas (semilla {semilla})")
