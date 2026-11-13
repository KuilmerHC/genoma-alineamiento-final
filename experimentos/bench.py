"""Harness común de medición (protocolo de la sección 4.3 del enunciado).

- Tiempo con time.perf_counter(), al menos 5 repeticiones por punto.
- Los datos se preparan ANTES de iniciar el cronómetro (fuera del tiempo medido).
- Se reporta mediana y dispersión (cuartiles Q1–Q3, mínimo y máximo).
- Pendiente log-log por mínimos cuadrados: si T(n) ≈ c·n^p entonces log T = log c + p·log n.
- El recolector de basura se apaga durante cada medición para reducir ruido.

Aquí se permiten numpy y sorted: es infraestructura de medición, no el núcleo de un componente.
"""

import csv
import gc
import os
import time

import matplotlib

matplotlib.use("Agg")  # sin ventana: solo genera archivos PNG
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

DIR_EXP = os.path.dirname(os.path.abspath(__file__))
DIR_RESULTADOS = os.path.join(DIR_EXP, "resultados")
DIR_GRAFICAS = os.path.join(DIR_EXP, "graficas")

# Paleta categórica validada (pares adyacentes distinguibles con daltonismo).
COLOR_MEDIDO = "#2a78d6"   # azul: datos medidos
COLOR_TEORIA = "#eb6834"   # naranja: predicción teórica
COLOR_EXTRA = "#1baf7a"    # aqua: tercera serie (p. ej. cota inferior)
TINTA = "#3b3b38"
GRILLA = "#e4e3dd"


def medir(funcion, preparar, repeticiones=5):
    """Ejecuta funcion(preparar()) `repeticiones` veces; devuelve la lista de tiempos (s).

    `preparar` genera una copia fresca de los datos en cada repetición, fuera del tiempo.
    """
    tiempos = []
    for _ in range(repeticiones):
        datos = preparar()
        gc.collect()
        gc.disable()
        try:
            t0 = time.perf_counter()
            funcion(datos)
            tiempos.append(time.perf_counter() - t0)
        finally:
            gc.enable()
    return tiempos


def resumen(tiempos):
    """Mediana y dispersión de una lista de tiempos."""
    t = np.asarray(tiempos)
    return {
        "mediana": float(np.median(t)),
        "q1": float(np.percentile(t, 25)),
        "q3": float(np.percentile(t, 75)),
        "min": float(t.min()),
        "max": float(t.max()),
    }


def pendiente_loglog(xs, ys):
    """Pendiente p del ajuste log(y) = p·log(x) + b (mínimos cuadrados)."""
    p, _ = np.polyfit(np.log(np.asarray(xs, float)), np.log(np.asarray(ys, float)), 1)
    return float(p)


def guardar_csv(nombre, filas):
    """Guarda una lista de diccionarios en experimentos/resultados/<nombre>.csv."""
    os.makedirs(DIR_RESULTADOS, exist_ok=True)
    ruta = os.path.join(DIR_RESULTADOS, nombre + ".csv")
    with open(ruta, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        w.writeheader()
        w.writerows(filas)
    return ruta


def _estilo(ax, titulo, xlabel, ylabel):
    ax.set_title(titulo, color=TINTA, fontsize=11, loc="left")
    ax.set_xlabel(xlabel, color=TINTA)
    ax.set_ylabel(ylabel, color=TINTA)
    ax.grid(True, which="major", color=GRILLA, linewidth=0.8)
    ax.tick_params(colors=TINTA, labelsize=9)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_color(GRILLA)


def grafica(nombre, xs, series, titulo, xlabel, ylabel, loglog=True, logx=False, nota=None,
            y_desde_cero=False):
    """Dibuja una o más series sobre UN solo eje y guarda graficas/<nombre>.png.

    series: lista de dicts con claves y, etiqueta, color y opcionales q1/q3 (barras de
    dispersión), estilo ('-' medido, '--' teoría) y marcador.
    y_desde_cero: para gráficas de razón (T/f(n)), el eje Y empieza en 0 para no exagerar
    variaciones pequeñas alrededor de la constante.
    """
    os.makedirs(DIR_GRAFICAS, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=150)
    for s in series:
        y = np.asarray(s["y"], float)
        if "q1" in s:
            err = [y - np.asarray(s["q1"]), np.asarray(s["q3"]) - y]
            ax.errorbar(xs, y, yerr=err, color=s["color"], linewidth=2, marker="o",
                        markersize=5, capsize=3, label=s["etiqueta"])
        else:
            ax.plot(xs, y, s.get("estilo", "-"), color=s["color"], linewidth=2,
                    marker=s.get("marcador"), markersize=5, label=s["etiqueta"])
    if loglog:
        ax.set_xscale("log")
        ax.set_yscale("log")
    elif logx:
        ax.set_xscale("log")
    if y_desde_cero:
        ax.set_ylim(bottom=0, top=max(max(s["y"]) for s in series) * 1.15)
    _estilo(ax, titulo, xlabel, ylabel)
    if len(series) > 1:
        ax.legend(frameon=False, fontsize=9, labelcolor=TINTA)
    if nota:
        ax.text(0.02, 0.97, nota, transform=ax.transAxes, va="top", fontsize=9, color=TINTA)
    fig.tight_layout()
    ruta = os.path.join(DIR_GRAFICAS, nombre + ".png")
    fig.savefig(ruta)
    plt.close(fig)
    return ruta


def tamanos_geometricos(base, cantidad):
    """[base, 2·base, 4·base, ...]: progresión geométrica (n se duplica)."""
    return [base * 2 ** i for i in range(cantidad)]
