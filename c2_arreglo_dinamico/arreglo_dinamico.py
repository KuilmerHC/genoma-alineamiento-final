"""C2: arreglo dinámico de posiciones de coincidencia (duplica su capacidad al llenarse).

Análisis amortizado por el método contable:
    Se cobran 3 créditos por append.
      - 1 crédito paga escribir el elemento.
      - 2 créditos quedan guardados en ese elemento.
    Al duplicar de capacidad c a 2c hay que copiar c elementos. Los c/2 elementos insertados
    desde la última duplicación guardan 2 créditos cada uno = c créditos, justo lo que cuesta
    la copia. El saldo nunca es negativo, así que el costo amortizado de append es O(1)
    y n appends cuestan O(n) (a lo sumo 3n).

    La clase lleva ese banco de créditos en `saldo` (y su mínimo histórico en `saldo_minimo`)
    para comprobar experimentalmente que nunca queda negativo.

Complejidad: agregar O(1) amortizado, Θ(n) en el peor caso (la operación que duplica);
acceso por índice Θ(1). Memoria: capacidad < 2n, es decir Θ(n).
"""

CREDITOS_POR_APPEND = 3


class ArregloDinamico:
    def __init__(self, capacidad=1):
        self.capacidad = capacidad
        self.datos = [None] * capacidad
        self.n = 0
        self.escrituras = 0   # escrituras de elementos nuevos
        self.copias = 0       # elementos copiados al redimensionar
        self.saldo = 0        # créditos cobrados − créditos gastados (método contable)
        self.saldo_minimo = 0

    def costo_total(self):
        return self.escrituras + self.copias

    def _redimensionar(self):
        nuevos = [None] * (2 * self.capacidad)
        for i in range(self.n):
            nuevos[i] = self.datos[i]
            self.copias += 1
            self.saldo -= 1   # cada copia se paga con 1 crédito guardado
        self.datos = nuevos
        self.capacidad *= 2

    def agregar(self, x):
        """Agrega x al final. O(1) amortizado; Θ(n) cuando hay que duplicar."""
        self.saldo += CREDITOS_POR_APPEND
        if self.n == self.capacidad:
            self._redimensionar()
        self.datos[self.n] = x
        self.escrituras += 1
        self.saldo -= 1       # la escritura del elemento nuevo cuesta 1 crédito
        self.n += 1
        self.saldo_minimo = min(self.saldo_minimo, self.saldo)

    def __getitem__(self, i):
        if i < 0 or i >= self.n:
            raise IndexError(i)
        return self.datos[i]

    def __len__(self):
        return self.n
