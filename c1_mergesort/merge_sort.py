"""C1: merge sort de k-mers (divide y vencerás) con contador de comparaciones.

Recurrencia (leída del código de `merge_sort`):
    T(n) = 2T(n/2) + Θ(n)    (dos llamadas sobre mitades + `mezclar` lineal)
Teorema Maestro: a=2, b=2, f(n)=n = Θ(n^{log_b a}) = Θ(n)  →  caso 2  →  T(n) = Θ(n log n).

Comparaciones:
    mejor caso  ≈ (n/2) log₂ n   (entrada ya ordenada: cada mezcla agota una mitad pronto)
    peor caso   ≈ n log₂ n − n + 1
    Toda la ejecución es Θ(n log n) en tiempo en los tres casos (las copias son lineales igual).
Cota inferior: cualquier ordenamiento por comparaciones necesita, en el peor caso,
    ⌈log₂ n!⌉ = Ω(n log n) comparaciones (altura mínima de un árbol de decisión con n! hojas).
Memoria extra: Θ(n) (listas auxiliares de la mezcla y slices).
Estable: con `<=` los elementos iguales conservan su orden relativo.
"""

comparaciones = 0  # contador global instrumentable (reiniciar antes de medir)


def mezclar(izq, der):
    """Mezcla dos listas ordenadas en una lista ordenada.

    Hace a lo sumo len(izq)+len(der)−1 comparaciones. Tiempo Θ(len(izq)+len(der)).
    """
    global comparaciones
    res = []
    i = j = 0
    while i < len(izq) and j < len(der):
        comparaciones += 1
        if izq[i] <= der[j]:
            res.append(izq[i])
            i += 1
        else:
            res.append(der[j])
            j += 1
    res.extend(izq[i:])
    res.extend(der[j:])
    return res


def merge_sort(a):
    """Devuelve una nueva lista con los elementos de `a` ordenados.

    Tiempo Θ(n log n) en mejor, peor y caso esperado. Memoria extra Θ(n).
    """
    if len(a) <= 1:
        return a
    mid = len(a) // 2
    return mezclar(merge_sort(a[:mid]), merge_sort(a[mid:]))


def entrada_peor_caso(ordenados):
    """Permutación de `ordenados` que obliga a merge sort a hacer n log₂ n − n + 1
    comparaciones (n potencia de 2): el peor caso.

    Idea: en cada mezcla, las dos mitades tienen que intercalarse hasta el final. Se logra
    "deshaciendo" una mezcla: los elementos en posiciones pares (en orden) van a la mitad
    izquierda y los de posiciones impares a la derecha, y se repite recursivamente en cada
    mitad. Ejemplo: [0..7] → [0, 4, 2, 6, 1, 5, 3, 7]. Tiempo Θ(n log n).
    """
    if len(ordenados) <= 1:
        return list(ordenados)
    pares = [ordenados[i] for i in range(0, len(ordenados), 2)]
    impares = [ordenados[i] for i in range(1, len(ordenados), 2)]
    return entrada_peor_caso(pares) + entrada_peor_caso(impares)


def reiniciar_contador():
    """Pone en cero el contador global de comparaciones."""
    global comparaciones
    comparaciones = 0
