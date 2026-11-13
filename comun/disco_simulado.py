"""Disco simulado (regla 4.1 del enunciado): el ÚNICO punto de acceso a datos en disco.

El archivo binario se ve como una secuencia de bloques de B registros. Cada registro es
(kmer: uint64, posicion: uint32) empaquetado con `struct` en 12 bytes (little-endian, sin
relleno). El bloque i ocupa los bytes [i·B·12, (i+1)·B·12). Solo el ÚLTIMO bloque del archivo
puede estar incompleto.

Cada `leer_bloque` y cada `escribir_bloque` suma 1 a su contador: eso es lo que se compara
con la fórmula de E/S del modelo de memoria externa. Ningún otro código abre el archivo.
"""

import os
import struct

REGISTRO = struct.Struct("<QI")  # kmer uint64 + posición uint32 = 12 bytes
TAM_REGISTRO = REGISTRO.size


class DiscoSimulado:
    def __init__(self, ruta, B, crear=True):
        """Abre (o crea vacío, si crear=True) el archivo `ruta` con bloques de B registros."""
        if B < 1:
            raise ValueError("B debe ser ≥ 1")
        self.ruta = ruta
        self.B = B
        self.bytes_por_bloque = B * TAM_REGISTRO
        self._f = open(ruta, "w+b" if crear else "r+b")
        self.lecturas = 0
        self.escrituras = 0

    # --- las dos operaciones contadas -------------------------------------------------
    def leer_bloque(self, i):
        """Devuelve la lista de registros (kmer, posicion) del bloque i. Cuenta 1 lectura."""
        self.lecturas += 1
        self._f.seek(i * self.bytes_por_bloque)
        datos = self._f.read(self.bytes_por_bloque)
        return list(REGISTRO.iter_unpack(datos))

    def escribir_bloque(self, i, registros):
        """Escribe hasta B registros en el bloque i. Cuenta 1 escritura."""
        if len(registros) > self.B:
            raise ValueError(f"un bloque admite a lo sumo {self.B} registros")
        self.escrituras += 1
        self._f.seek(i * self.bytes_por_bloque)
        self._f.write(b"".join(REGISTRO.pack(k, p) for k, p in registros))

    # --- metadatos (no son transferencias de bloques, no se cuentan) --------------------
    def num_registros(self):
        self._f.seek(0, os.SEEK_END)
        return self._f.tell() // TAM_REGISTRO

    def num_bloques(self):
        return -(-self.num_registros() // self.B)  # techo de N/B

    @property
    def total_es(self):
        return self.lecturas + self.escrituras

    def reiniciar_contadores(self):
        self.lecturas = 0
        self.escrituras = 0

    def cerrar(self):
        self._f.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.cerrar()


# --- utilidades para preparar y verificar datos (fuera de la medición) ---------------------
def escribir_registros(disco, registros):
    """Vuelca una lista de registros al disco, bloque por bloque (usa escribir_bloque)."""
    B = disco.B
    for b, inicio in enumerate(range(0, len(registros), B)):
        disco.escribir_bloque(b, registros[inicio:inicio + B])


def leer_todo(disco):
    """Lee todos los registros del disco, bloque por bloque (usa leer_bloque)."""
    res = []
    for b in range(disco.num_bloques()):
        res.extend(disco.leer_bloque(b))
    return res
