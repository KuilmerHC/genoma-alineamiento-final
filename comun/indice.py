"""Índice de semillas: une C1 (merge sort) y C2 (arreglo dinámico de posiciones).

Construcción, para un genoma de largo n y u k-mers distintos:
  1. extraer los n−k+1 pares (kmer_codificado, posicion)          Θ(n)
  2. ordenarlos con merge sort (C1)                               Θ(n log n)
  3. recorrer la lista ordenada: cada k-mer distinto se guarda una vez en `claves` y sus
     posiciones se agregan a un ArregloDinamico (C2)              Θ(n) amortizado
  Total Θ(n log n), dominado por el ordenamiento.
Consulta: búsqueda binaria sobre `claves` (ordenadas)            Θ(log u) ⊆ O(log n)
"""

from c1_mergesort.merge_sort import merge_sort
from c2_arreglo_dinamico.arreglo_dinamico import ArregloDinamico
from comun.kmers import kmers_codificados


class IndiceSemillas:
    def __init__(self, genoma, k, ordenar=merge_sort):
        """Construye el índice. `ordenar` permite cambiar C1 por C4 (quicksort) después."""
        self.k = k
        self.claves = []       # k-mers distintos, en orden creciente
        self.posiciones = []   # posiciones[t] = ArregloDinamico con las posiciones de claves[t]
        for kmer, pos in ordenar(kmers_codificados(genoma, k)):
            if not self.claves or self.claves[-1] != kmer:
                self.claves.append(kmer)
                self.posiciones.append(ArregloDinamico())
            self.posiciones[-1].agregar(pos)

    def buscar(self, kmer):
        """Posiciones del k-mer (codificado) en el genoma, o una lista vacía.

        Búsqueda binaria: invariante "si el k-mer está, su índice está en [lo, hi]".
        Cada iteración parte el intervalo a la mitad → T(u) = T(u/2) + Θ(1) = Θ(log u).
        """
        lo, hi = 0, len(self.claves) - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            if self.claves[mid] == kmer:
                return self.posiciones[mid]
            if self.claves[mid] < kmer:
                lo = mid + 1
            else:
                hi = mid - 1
        return []
