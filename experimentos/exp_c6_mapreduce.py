"""Experimento C6: MapReduce con multiprocessing. Speedup frente a la ley de Amdahl.

Hipótesis H6:
  (a) el trabajo crece linealmente con el número de lecturas R (pendiente log-log ≈ 1, p = 1);
  (b) el speedup medido S(p) = T₁/T_p sigue a Amdahl, S(p) = 1/(f + (1 − f)/p), usando la
      fracción serial f MEDIDA (no ajustada):
        escenario "solo alineamiento": índice ya construido; f = (arranque + reduce)/T₁, pequeña;
        escenario "sistema completo":  incluye construir el índice (serial); f grande, y el
                                       speedup queda acotado por 1/f.
  (c) más allá de 8 procesos no se gana: la máquina tiene 8 núcleos físicos (16 hilos SMT).

Diagnóstico: métrica de Karp-Flatt, la fracción serial EXPERIMENTAL que implica cada
speedup medido, e(p) = (1/S − 1/p) / (1 − 1/p). Si e(p) es casi constante, el sobrecosto se
comporta como una parte serial fija (Amdahl con f = e); si crece con p, es un sobrecosto
que aumenta con el número de procesos (comunicación, contención).

Uso: python -m experimentos.exp_c6_mapreduce
"""

import os

from c6_mapreduce.mapreduce import alinear_mapreduce
from comun.indice import IndiceSemillas
from datos.generar import generar_genoma, generar_lecturas
from experimentos.bench import (COLOR_EXTRA, COLOR_MEDIDO, COLOR_TEORIA, grafica, guardar_csv,
                                medir, pendiente_loglog, resumen)

LARGO_GENOMA = 10**6
R = 16_000
PROCESOS = (1, 2, 4, 8, 16)


def amdahl(f, p):
    return 1.0 / (f + (1.0 - f) / p)


def main(repeticiones=5):
    genoma = generar_genoma(LARGO_GENOMA, semilla=42)
    lecturas = [l for _, l in generar_lecturas(genoma, R, semilla=42, tasa_indel=0.01)]
    print(f"núcleos lógicos: {os.cpu_count()}")

    # Fase serial 1: construir el índice (se mide aparte).
    r_indice = resumen(medir(lambda g: IndiceSemillas(g, 12), lambda: genoma, repeticiones))
    indice = IndiceSemillas(genoma, 12)
    print(f"índice: {r_indice['mediana']:.2f}s")

    # (a) Trabajo lineal en R con p = 1.
    filas_r = []
    for r in (1000, 2000, 4000, 8000, 16000):
        res = resumen(medir(lambda L: alinear_mapreduce(genoma, L, 1, indice=indice),
                            lambda: lecturas[:r], repeticiones))
        filas_r.append({"R": r, **res, "ms_por_lectura": res["mediana"] / r * 1e3})
        print(f"R={r:>6}  T1={res['mediana']:.2f}s  {filas_r[-1]['ms_por_lectura']:.3f} ms/lectura")
    p_r = pendiente_loglog([f["R"] for f in filas_r], [f["mediana"] for f in filas_r])
    print(f"pendiente log-log T1 vs R = {p_r:.3f}")
    guardar_csv("c6_trabajo_vs_R", filas_r)

    # (b) Speedup con R fijo.
    filas_p = []
    for p in PROCESOS:
        fases = []

        def correr(L, p=p):
            _, _, t = alinear_mapreduce(genoma, L, p, indice=indice)
            fases.append(t)

        res = resumen(medir(correr, lambda: lecturas, repeticiones))
        mediana = lambda clave: sorted(t[clave] for t in fases)[len(fases) // 2]  # noqa: E731
        filas_p.append({"p": p, **res, "arranque": mediana("arranque"), "map": mediana("map"),
                        "reduce": mediana("reduce")})
        print(f"p={p:>2}  T={res['mediana']:.2f}s  map={filas_p[-1]['map']:.2f}s "
              f"arranque={filas_p[-1]['arranque']:.3f}s")

    T1 = filas_p[0]["mediana"]
    T_indice = r_indice["mediana"]
    # Fracción serial medida. Con p = 1 no hay arranque de procesos; se toma el arranque
    # medido con p = 2 como el costo serial de crear y cerrar procesos.
    serial_sin_indice = filas_p[1]["arranque"] + filas_p[0]["reduce"]
    f_alin = serial_sin_indice / T1
    f_total = (T_indice + serial_sin_indice) / (T_indice + T1)
    for fila in filas_p:
        p = fila["p"]
        fila["speedup_alineamiento"] = T1 / fila["mediana"]
        fila["amdahl_alineamiento"] = amdahl(f_alin, p)
        fila["speedup_completo"] = (T_indice + T1) / (T_indice + fila["mediana"])
        fila["amdahl_completo"] = amdahl(f_total, p)
        fila["eficiencia_alineamiento"] = fila["speedup_alineamiento"] / p
        S = fila["speedup_alineamiento"]
        fila["karp_flatt"] = (1 / S - 1 / p) / (1 - 1 / p) if p > 1 else float("nan")
    # Fracción serial efectiva: promedio de Karp-Flatt en p = 2, 4, 8 (núcleos físicos).
    kf = [f["karp_flatt"] for f in filas_p if 1 < f["p"] <= 8]
    f_ef = sum(kf) / len(kf)
    for fila in filas_p:
        fila["amdahl_f_efectiva"] = amdahl(f_ef, fila["p"])
    print(f"f (solo alineamiento) = {f_alin:.4f}  →  S(∞) ≤ {1 / f_alin:.1f}")
    print(f"f (sistema completo)  = {f_total:.4f}  →  S(∞) ≤ {1 / f_total:.2f}")
    print(f"f efectiva (Karp-Flatt, p = 2..8) = {f_ef:.4f}")
    for f in filas_p:
        print(f"p={f['p']:>2}  S_alin={f['speedup_alineamiento']:.2f} (Amdahl {f['amdahl_alineamiento']:.2f})"
              f"  S_completo={f['speedup_completo']:.2f} (Amdahl {f['amdahl_completo']:.2f})"
              f"  eficiencia={f['eficiencia_alineamiento']:.2f}  Karp-Flatt={f['karp_flatt']:.4f}"
              f"  Amdahl(f_ef)={f['amdahl_f_efectiva']:.2f}")
    for f in filas_p:
        f["f_alineamiento"] = f_alin
        f["f_completo"] = f_total
        f["f_efectiva"] = f_ef
        f["T_indice"] = T_indice
    guardar_csv("c6_speedup", filas_p)

    ps = list(PROCESOS)
    grafica("c6_speedup_alineamiento", ps, [
        {"y": [f["speedup_alineamiento"] for f in filas_p], "etiqueta": "medido",
         "color": COLOR_MEDIDO, "marcador": "o"},
        {"y": [f["amdahl_alineamiento"] for f in filas_p],
         "etiqueta": f"Amdahl, f = {f_alin:.3f} (serial medida)", "color": COLOR_TEORIA,
         "estilo": "--"},
        {"y": [f["amdahl_f_efectiva"] for f in filas_p],
         "etiqueta": f"Amdahl, f = {f_ef:.3f} (Karp-Flatt)", "color": COLOR_EXTRA, "estilo": ":"},
    ], "C6: speedup del alineamiento (índice ya construido)", "procesos p", "speedup T₁/T_p",
        loglog=False, y_desde_cero=True, nota="8 núcleos físicos, 16 hilos")
    grafica("c6_speedup_completo", ps, [
        {"y": [f["speedup_completo"] for f in filas_p], "etiqueta": "medido",
         "color": COLOR_MEDIDO, "marcador": "o"},
        {"y": [f["amdahl_completo"] for f in filas_p],
         "etiqueta": f"Amdahl, f = {f_total:.3f} (medida)", "color": COLOR_TEORIA, "estilo": "--"},
        {"y": [1 / f_total] * len(ps), "etiqueta": f"cota 1/f = {1 / f_total:.2f}",
         "color": COLOR_EXTRA, "estilo": ":"},
    ], "C6: speedup del sistema completo (con índice serial)", "procesos p", "speedup",
        loglog=False, y_desde_cero=True)
    grafica("c6_trabajo_vs_R", [f["R"] for f in filas_r], [
        {"y": [f["mediana"] for f in filas_r], "q1": [f["q1"] for f in filas_r],
         "q3": [f["q3"] for f in filas_r], "etiqueta": "T₁ medido", "color": COLOR_MEDIDO},
    ], "C6: trabajo con p = 1 según el número de lecturas", "R (lecturas)", "tiempo (s)",
        nota=f"pendiente = {p_r:.3f}")
    return filas_r, filas_p


if __name__ == "__main__":
    main()
