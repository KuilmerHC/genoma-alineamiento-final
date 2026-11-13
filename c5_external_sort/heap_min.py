"""Min-heap binario implementado por nosotros (el enunciado prohíbe heapq en el núcleo).

Se guarda en una lista como árbol binario completo: los hijos del nodo i están en 2i+1 y
2i+2, y su padre en (i−1)//2. Propiedad de heap: cada nodo es ≤ que sus hijos, así que el
mínimo está en la raíz (posición 0).

    insertar:     O(log k)  (se agrega al final y "sube" mientras sea menor que su padre)
    extraer_min:  O(log k)  (el último pasa a la raíz y "baja" hacia el hijo menor)
    con k = número de elementos en el heap.
"""


class MinHeap:
    def __init__(self):
        self._a = []

    def __len__(self):
        return len(self._a)

    def insertar(self, x):
        a = self._a
        a.append(x)
        i = len(a) - 1
        while i > 0:
            padre = (i - 1) // 2
            if a[i] < a[padre]:
                a[i], a[padre] = a[padre], a[i]
                i = padre
            else:
                break

    def extraer_min(self):
        a = self._a
        minimo = a[0]
        ultimo = a.pop()
        if a:
            a[0] = ultimo
            i, n = 0, len(a)
            while True:
                izq, der = 2 * i + 1, 2 * i + 2
                menor = i
                if izq < n and a[izq] < a[menor]:
                    menor = izq
                if der < n and a[der] < a[menor]:
                    menor = der
                if menor == i:
                    break
                a[i], a[menor] = a[menor], a[i]
                i = menor
        return minimo
