"""C4: quicksort aleatorizado con partición de 3 vías (Dijkstra), con contador de comparaciones.

Idea: se elige un pivote AL AZAR y se reordena el tramo en tres zonas
    [ < pivote | = pivote | > pivote ]
Luego se ordenan recursivamente solo las zonas < y >. La zona = ya quedó en su lugar final,
lo que importa mucho con k-mers repetidos.

Qué se cuenta: cada comparación de un elemento contra el pivote cuenta UNA vez (es una
comparación de tres resultados: <, = o >). Es lo mismo que cuenta el análisis clásico.

Caso esperado (claves distintas). Sean z_1 < z_2 < … < z_n los elementos en orden. z_i y z_j
se comparan si y solo si el primer pivote elegido entre {z_i, …, z_j} es z_i o z_j: eso
ocurre con probabilidad 2/(j−i+1). Por linealidad de la esperanza:
    E[comparaciones] = Σ_{i<j} 2/(j−i+1) = 2(n+1)·H_n − 4n ≈ 2n ln n ≈ 1,39·n·log₂ n
    (H_n = 1 + 1/2 + … + 1/n ≈ ln n + 0,577).
La esperanza es sobre las elecciones aleatorias del algoritmo, NO sobre la entrada: vale
para cualquier entrada, incluida una ya ordenada.

Peor caso: Θ(n²), si cada pivote resulta ser el mínimo o el máximo. Con pivote aleatorio
la probabilidad es despreciable (≤ 2^(n−1)/n! para la secuencia exacta de pivotes malos).
Mejor caso: Θ(n log n) con claves distintas (pivote siempre la mediana); Θ(n) si todas las
claves son iguales (una sola partición deja todo en la zona =).
Con u claves distintas, esperado O(n log u).
Memoria extra: O(log n) en el peor caso, porque se hace recursión solo sobre la zona más
pequeña (tamaño ≤ la mitad) y la más grande se procesa con un ciclo.
No es estable.
"""

import random

comparaciones = 0  # contador global instrumentable (reiniciar antes de medir)


def reiniciar_contador():
    """Pone en cero el contador global de comparaciones."""
    global comparaciones
    comparaciones = 0


def particion_3_vias(a, lo, hi, rng):
    """Particiona a[lo..hi] (inclusive) alrededor de un pivote aleatorio.

    Devuelve (lt, gt) tales que
        a[lo..lt-1] < pivote,  a[lt..gt] == pivote,  a[gt+1..hi] > pivote.
    Invariante del ciclo: a[lo..lt-1] < p, a[lt..i-1] == p, a[gt+1..hi] > p, y a[i..gt]
    está sin revisar. Cada paso revisa a[i] (una comparación contra p) y achica la zona sin
    revisar en 1. Tiempo Θ(hi − lo + 1).
    """
    global comparaciones
    k = rng.randint(lo, hi)
    a[lo], a[k] = a[k], a[lo]  # el pivote aleatorio se lleva al inicio
    p = a[lo]
    lt, i, gt = lo, lo + 1, hi
    while i <= gt:
        comparaciones += 1  # comparar a[i] contra p (resultado <, = o >)
        x = a[i]
        if x < p:
            a[lt], a[i] = x, a[lt]
            lt += 1
            i += 1
        elif x > p:
            a[gt], a[i] = x, a[gt]
            gt -= 1
        else:
            i += 1
    return lt, gt


def quicksort(a, semilla=None):
    """Ordena la lista `a` EN SU LUGAR y la devuelve.

    `semilla` fija los pivotes aleatorios (para experimentos reproducibles).
    Tiempo esperado Θ(n log n) para cualquier entrada; peor caso Θ(n²) con prob. despreciable.
    """
    rng = random.Random(semilla)
    _quicksort(a, 0, len(a) - 1, rng)
    return a


def _quicksort(a, lo, hi, rng):
    while lo < hi:
        lt, gt = particion_3_vias(a, lo, hi, rng)
        # Recursión sobre la zona más pequeña; la grande se sigue en el mismo ciclo.
        if lt - lo < hi - gt:
            _quicksort(a, lo, lt - 1, rng)
            lo = gt + 1
        else:
            _quicksort(a, gt + 1, hi, rng)
            hi = lt - 1


def comparaciones_esperadas(n):
    """2(n+1)·H_n − 4n: comparaciones esperadas con n claves distintas."""
    if n < 2:
        return 0.0
    h = 0.0
    for i in range(1, n + 1):
        h += 1.0 / i
    return 2 * (n + 1) * h - 4 * n
