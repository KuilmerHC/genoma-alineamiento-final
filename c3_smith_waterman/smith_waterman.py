"""C3: Smith-Waterman, alineamiento local por programación dinámica.

Problema: dadas una lectura `a` (largo m) y una región del genoma `b` (largo n), encontrar
el par de subcadenas a[i'..i) y b[j'..j) cuyo alineamiento tiene el puntaje máximo.

Puntajes (sistema lineal de huecos):
    coincidencia (misma base)       +COINCIDENCIA   (por defecto +2)
    sustitución (bases distintas)   DIFERENCIA      (por defecto −1)
    hueco (inserción o borrado)     HUECO           (por defecto −2)

Recurrencia. H[i][j] = mejor puntaje de un alineamiento local que TERMINA en a[i-1] y b[j-1]:
    H[i][0] = H[0][j] = 0
    H[i][j] = max( 0,                                   empezar de nuevo (alineamiento vacío)
                   H[i-1][j-1] + s(a[i-1], b[j-1]),     alinear las dos bases
                   H[i-1][j]   + HUECO,                 a[i-1] contra un hueco
                   H[i][j-1]   + HUECO )                b[j-1] contra un hueco
    Respuesta: max_{i,j} H[i][j]. El alineamiento se reconstruye retrocediendo desde esa
    celda hasta llegar a un 0.

Subestructura óptima: si un alineamiento óptimo que termina en (i, j) usa su última columna
(diagonal, arriba o izquierda), lo que queda antes es un alineamiento óptimo que termina en
(i-1, j-1), (i-1, j) o (i, j-1). Si no lo fuera, se reemplazaría por uno mejor y el total
mejoraría, lo cual es una contradicción (argumento de cortar y pegar).
Subproblemas superpuestos: H[i-1][j-1] lo usan H[i][j], H[i-1][j] y H[i][j-1]. Una recursión
ingenua lo recalcularía un número exponencial de veces; la tabla lo calcula una sola vez.

Complejidad: tiempo Θ(mn) en mejor, peor y caso esperado (se llenan las m·n celdas, cada una
en O(1)); retroceso O(m + n). Memoria Θ(mn) por la tabla completa.
"""

from typing import NamedTuple

COINCIDENCIA = 2
DIFERENCIA = -1
HUECO = -2


class Alineamiento(NamedTuple):
    puntaje: int
    inicio_a: int   # subcadena alineada de a: a[inicio_a:fin_a]
    fin_a: int
    inicio_b: int   # subcadena alineada de b: b[inicio_b:fin_b]
    fin_b: int
    alin_a: str     # a alineada, con '-' en los huecos
    alin_b: str     # b alineada, con '-' en los huecos
    celdas: int     # celdas de la tabla calculadas (= m·n), para validar Θ(mn)


def llenar_tabla(a, b, coincidencia=COINCIDENCIA, diferencia=DIFERENCIA, hueco=HUECO):
    """Llena la tabla H de (m+1)×(n+1). Devuelve (H, puntaje_max, i_max, j_max). Θ(mn)."""
    m, n = len(a), len(b)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    mejor, i_max, j_max = 0, 0, 0
    for i in range(1, m + 1):
        fila, fila_ant = H[i], H[i - 1]
        ai = a[i - 1]
        for j in range(1, n + 1):
            s = coincidencia if ai == b[j - 1] else diferencia
            h = max(0,
                    fila_ant[j - 1] + s,   # diagonal
                    fila_ant[j] + hueco,   # arriba
                    fila[j - 1] + hueco)   # izquierda
            fila[j] = h
            if h > mejor:
                mejor, i_max, j_max = h, i, j
    return H, mejor, i_max, j_max


def smith_waterman(a, b, coincidencia=COINCIDENCIA, diferencia=DIFERENCIA, hueco=HUECO):
    """Alineamiento local óptimo de `a` contra `b`. Tiempo y memoria Θ(len(a)·len(b))."""
    H, mejor, i, j = llenar_tabla(a, b, coincidencia, diferencia, hueco)
    fin_a, fin_b = i, j
    alin_a, alin_b = [], []
    # Retroceso: en cada paso se elige la celda de la que vino el valor H[i][j].
    # Se para al llegar a 0, porque ahí empieza el alineamiento local.
    while i > 0 and j > 0 and H[i][j] > 0:
        s = coincidencia if a[i - 1] == b[j - 1] else diferencia
        if H[i][j] == H[i - 1][j - 1] + s:
            alin_a.append(a[i - 1])
            alin_b.append(b[j - 1])
            i, j = i - 1, j - 1
        elif H[i][j] == H[i - 1][j] + hueco:
            alin_a.append(a[i - 1])
            alin_b.append("-")
            i -= 1
        else:  # necesariamente H[i][j] == H[i][j-1] + hueco
            alin_a.append("-")
            alin_b.append(b[j - 1])
            j -= 1
    return Alineamiento(
        puntaje=mejor,
        inicio_a=i, fin_a=fin_a,
        inicio_b=j, fin_b=fin_b,
        alin_a="".join(reversed(alin_a)),
        alin_b="".join(reversed(alin_b)),
        celdas=len(a) * len(b),
    )


def puntaje_dos_filas(a, b, coincidencia=COINCIDENCIA, diferencia=DIFERENCIA, hueco=HUECO):
    """Solo el puntaje máximo, guardando dos filas de la tabla en vez de toda la tabla.

    H[i][j] depende únicamente de la fila i−1 y de la fila i, así que basta con la fila
    anterior y la actual. Tiempo Θ(mn), igual que antes; memoria Θ(n) en vez de Θ(mn).
    Precio: sin la tabla completa no se puede hacer el retroceso (no hay alineamiento).
    """
    n = len(b)
    anterior = [0] * (n + 1)
    mejor = 0
    for i in range(1, len(a) + 1):
        actual = [0] * (n + 1)
        ai = a[i - 1]
        for j in range(1, n + 1):
            s = coincidencia if ai == b[j - 1] else diferencia
            h = max(0, anterior[j - 1] + s, anterior[j] + hueco, actual[j - 1] + hueco)
            actual[j] = h
            if h > mejor:
                mejor = h
        anterior = actual
    return mejor


def puntaje_de_alineamiento(alin_a, alin_b, coincidencia=COINCIDENCIA,
                            diferencia=DIFERENCIA, hueco=HUECO):
    """Puntaje de un alineamiento ya armado (columna por columna). Se usa para verificar."""
    total = 0
    for x, y in zip(alin_a, alin_b):
        if x == "-" or y == "-":
            total += hueco
        elif x == y:
            total += coincidencia
        else:
            total += diferencia
    return total
