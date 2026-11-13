"""Codificación de k-mers en 2 bits por base.

Cada base se codifica como A=0, C=1, G=2, T=3, y un k-mer s[0..k-1] como el entero
    x = Σ cod(s[i]) · 4^(k-1-i)
es decir, la base de más a la izquierda queda en los bits más significativos.
Consecuencias:
  - Un k-mer con k ≤ 31 cabe en 62 bits (entra en uint64, que es lo que usa C5 en disco).
  - Comparar dos k-mers es comparar dos enteros: O(1) en vez de O(k) con strings.
  - El orden de los enteros coincide con el orden alfabético de los strings (A<C<G<T),
    porque es un número en base 4 con los dígitos en el mismo orden.
"""

COD = {"A": 0, "C": 1, "G": 2, "T": 3}
BASES = "ACGT"
K_MAX = 31


def codificar_kmer(s):
    """Codifica el string `s` (solo A, C, G, T) como entero. Tiempo Θ(k)."""
    if len(s) > K_MAX:
        raise ValueError(f"k={len(s)} > {K_MAX}: no cabe en 64 bits")
    x = 0
    for base in s:
        x = (x << 2) | COD[base]  # correr 2 bits a la izquierda y poner la base al final
    return x


def decodificar_kmer(x, k):
    """Inverso de `codificar_kmer`: entero → string de largo k. Tiempo Θ(k)."""
    letras = []
    for _ in range(k):
        letras.append(BASES[x & 3])  # los 2 bits más bajos son la última base
        x >>= 2
    return "".join(reversed(letras))


def kmers_codificados(texto, k):
    """Lista de (kmer_codificado, posicion) para todas las ventanas de largo k del texto.

    Ventana deslizante: el k-mer en i+1 se obtiene del de i corriendo 2 bits, agregando
    la base nueva y borrando la base que sale con una máscara de 2k bits.
    Tiempo Θ(n) para un texto de largo n (en vez de Θ(nk) codificando cada ventana).
    """
    if not 1 <= k <= K_MAX:
        raise ValueError(f"k debe estar entre 1 y {K_MAX}")
    n = len(texto)
    if n < k:
        return []
    mascara = (1 << (2 * k)) - 1
    x = codificar_kmer(texto[:k])
    res = [(x, 0)]
    for i in range(1, n - k + 1):
        x = ((x << 2) | COD[texto[i + k - 1]]) & mascara
        res.append((x, i))
    return res
