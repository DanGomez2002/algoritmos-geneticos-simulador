"""
graficas.py
-----------
Gráficas del simulador. Cada función dibuja sobre unos ejes que recibe, así
se usan igual en la ventana interactiva y en generar_resultados.py.
"""

import textwrap

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Arc, Circle, Rectangle
from matplotlib.textpath import TextPath
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Affine2D

import plantilla as pl

VERDE_CANCHA = "#3f8f4f"
VERDE_OSCURO = "#367c44"
ROJO = "#c53030"
AZUL = "#2b6cb0"
GRIS = "#718096"
VERDE = "#2f855a"

# Un color fijo por jugador (para seguirlo en la tabla de genes).
COLORES_JUGADOR = [plt.get_cmap("tab20")(i) for i in range(20)][:pl.N_JUGADORES]
# Color del borde según el tipo de rol del puesto.
COLOR_ROL = {"POR": "#ecc94b", "DFC": "#3182ce", "LAT": "#3182ce", "CAR": "#3182ce",
             "MCD": "#38a169", "MC": "#38a169", "MED": "#38a169", "EXT": "#e53e3e", "DC": "#e53e3e"}

# Rol natural de cada jugador expresado como grupo, para marcar quién juega fuera de posición.
_GRUPO = {"POR": "POR", "DFC": "DEF", "LD": "DEF", "LI": "DEF", "MCD": "MED", "MC": "MED",
          "MCO": "MED", "ED": "ATA", "EI": "ATA", "DC": "ATA"}
_GRUPO_ROL = {"POR": "POR", "DFC": "DEF", "LAT": "DEF", "CAR": "DEF", "MCD": "MED", "MC": "MED",
              "MED": "MED", "EXT": "ATA", "DC": "ATA"}


def fuera_de_posicion(jugador, rol):
    """True si el jugador juega fuera de su línea natural.

    Un extremo (ED/EI) de volante por banda no cuenta como fuera de posición:
    en un 4-4-2 ese es su puesto habitual.
    """
    if pl.NATURAL[jugador] in ("ED", "EI") and rol == "MED":
        return False
    return _GRUPO[pl.NATURAL[jugador]] != _GRUPO_ROL[rol]


# ---------------------------------------------------------------------------
# Cancha
# ---------------------------------------------------------------------------

def _lineas_cancha(ax):
    ax.add_patch(Rectangle((-4, -4), 108, 108, color=VERDE_CANCHA, zorder=0))
    for i in range(0, 100, 10):                         # franjas del césped
        if (i // 10) % 2:
            ax.add_patch(Rectangle((0, i), 100, 10, color=VERDE_OSCURO, zorder=0))
    kw = dict(color="white", lw=1.2, zorder=1)
    ax.plot([0, 100, 100, 0, 0], [0, 0, 100, 100, 0], **kw)
    ax.plot([0, 100], [50, 50], **kw)
    ax.add_patch(Circle((50, 50), 9, fill=False, ec="white", lw=1.2, zorder=1))
    for y0, s in ((0, 1), (100, -1)):
        ax.plot([20, 20, 80, 80], [y0, y0 + s * 16, y0 + s * 16, y0], **kw)
        ax.plot([37, 37, 63, 63], [y0, y0 + s * 6, y0 + s * 6, y0], **kw)
        ax.add_patch(Arc((50, y0 + s * 11), 18, 12, theta1=35 if s == 1 else 215,
                         theta2=145 if s == 1 else 325, color="white", lw=1.2, zorder=1))


def dibujar_cancha(ax, problema, cromosoma, titulo=None):
    """La alineación del cromosoma sobre la cancha, con las líneas de química."""
    ax.clear()
    _lineas_cancha(ax)
    c = np.asarray(cromosoma)
    puestos = problema.puestos
    # Química entre vecinos (verde = suma, rojo = resta).
    for i, k in problema.pares:
        q = pl.QUIMICA[c[i], c[k]]
        if q != 0:
            (x1, y1), (x2, y2) = puestos[i][3:5], puestos[k][3:5]
            ax.plot([x1, x2], [y1, y2], color="#9ae6b4" if q > 0 else "#feb2b2",
                    lw=1.5 + 1.5 * abs(q), alpha=0.95, zorder=2, solid_capstyle="round")
    lesionados = set(problema.lesionados.tolist())
    for p, (sigla, rol, lado, x, y) in enumerate(puestos):
        j = int(c[p])
        fuera = fuera_de_posicion(j, rol)
        ax.add_patch(Circle((x, y), 4.6, fc="white", ec=COLOR_ROL[rol], lw=3, zorder=3))
        if fuera:
            ax.add_patch(Circle((x, y), 6.0, fill=False, ec="#f6ad55", lw=1.6, ls="--", zorder=3))
        ax.text(x, y + 0.3, str(j + 1), ha="center", va="center", fontsize=8.5, weight="bold",
                color=ROJO if j in lesionados else "#1a202c", zorder=4)
        ax.text(x, y - 7.6, pl.CORTOS[j], ha="center", va="center", fontsize=7.2, color="white",
                weight="bold", zorder=4)
        ax.text(x, y + 6.6, f"{sigla} {problema.R[j, p]:.0f}", ha="center", va="center",
                fontsize=6.2, color="#e2e8f0", zorder=4)
    banco = ", ".join(f"{j + 1} {pl.CORTOS[j]}" for j in c[pl.N_TITULARES:])
    ax.text(50, -6, textwrap.fill("Banco: " + banco, 60), ha="center", va="top", fontsize=6.8)
    ax.set_xlim(-4, 104)
    ax.set_ylim(-14, 104)
    ax.set_aspect("auto")
    ax.axis("off")
    ax.set_title(titulo or "Mejor alineación", fontsize=9.5)


# ---------------------------------------------------------------------------
# Convergencia y diversidad
# ---------------------------------------------------------------------------

def dibujar_convergencia(ax, historial, referencia=None, gen_marcada=None, ancho_x=None):
    """Mejor puntaje y promedio de la población por generación."""
    ax.clear()
    g = np.arange(len(historial.mejor))
    ax.plot(g, historial.media, color=GRIS, lw=1.2, label="Promedio de la población")
    ax.plot(g, historial.mejor, color=AZUL, lw=1, alpha=0.6, label="Mejor de la generación")
    ax.plot(g, historial.mejor_global, color=ROJO, lw=2, label="Mejor encontrado (élite)")
    if referencia is not None:
        ax.axhline(referencia, color=VERDE, ls="--", lw=1.2,
                   label="Asignación húngara (ignora la química)")
    if gen_marcada is not None:
        ax.axvline(gen_marcada, color="k", lw=0.8, ls=":")
    valores = np.concatenate([historial.mejor_global, historial.media[-max(1, len(g) // 2):]])
    bajo = min(np.percentile(historial.media, 20), np.min(valores))
    alto = max(np.max(historial.mejor_global), referencia or -np.inf)
    ax.set_ylim(bajo - 0.04 * (alto - bajo), alto + 0.12 * (alto - bajo + 1))
    if ancho_x:
        ax.set_xlim(0, ancho_x)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("Generación", fontsize=8)
    ax.set_ylabel("Puntaje del equipo", fontsize=8)
    ax.set_title("Evolución del puntaje", fontsize=9.5)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6.5, loc="lower right")


# ---------------------------------------------------------------------------
# Genes de la población
# ---------------------------------------------------------------------------

def _numeros_en_celdas(ax, P, fontsize):
    """Escribe el número de cada jugador en su celda de la tabla de genes.

    Con 30 cromosomas x 18 genes serían 540 textos, y Matplotlib tarda cerca
    de 1 ms en dibujar cada uno (medio segundo por generación). En su lugar,
    cada número (1 a 18) se convierte UNA vez en una figura vectorial y se
    estampa en todas las celdas donde aparece: 18 objetos en vez de 540, con
    el mismo aspecto y unas diez veces más rápido.
    """
    # puntos -> pulgadas -> píxeles (sigue la resolución de la figura, también al guardar)
    escala = Affine2D().scale(1 / 72) + ax.figure.dpi_scale_trans
    for j in range(pl.N_JUGADORES):
        filas, cols = np.nonzero(P == j)
        if filas.size == 0:
            continue
        forma = TextPath((0, 0), str(j + 1), prop=FontProperties(size=fontsize))
        caja = forma.get_extents()
        forma = forma.transformed(Affine2D().translate(-(caja.x0 + caja.x1) / 2, -(caja.y0 + caja.y1) / 2))
        ax.add_collection(PathCollection([forma], offsets=np.column_stack([cols, filas]),
                                         offset_transform=ax.transData, transform=escala,
                                         facecolors="k", edgecolors="none"))


def dibujar_genes(ax, problema, poblacion, puntajes, n_elite=0, max_filas=15, generacion=None):
    """Cada fila es un cromosoma (del mejor al peor) y cada columna un gen.

    El color de la celda es el jugador que ocupa ese gen. Al principio la
    tabla es un mosaico de colores; con las generaciones las columnas se
    vuelven de un solo color: la población converge a una alineación.
    Se muestran solo los `max_filas` mejores para que los números se lean
    también en pantallas pequeñas.
    """
    ax.clear()
    total = len(poblacion)
    P = np.asarray(poblacion)[:max_filas]
    filas, L = P.shape
    imagen = np.array([[COLORES_JUGADOR[j] for j in fila] for fila in P])
    ax.imshow(imagen, aspect="auto", interpolation="nearest")
    _numeros_en_celdas(ax, P, fontsize=6)
    ax.axvline(pl.N_TITULARES - 0.5, color="k", lw=1.5)
    if n_elite:
        ax.add_patch(Rectangle((-0.5, -0.5), L, min(n_elite, filas), fill=False, ec=ROJO, lw=1.8))
    siglas = [p[0] for p in problema.puestos] + ["B"] * (L - pl.N_TITULARES)
    ax.set_xticks(range(L))
    ax.set_xticklabels(siglas, fontsize=6)
    ax.set_yticks([0, filas - 1])
    ax.set_yticklabels([f"1.º {puntajes[0]:.0f}", f"{filas}.º {puntajes[filas - 1]:.0f}"], fontsize=6.5)
    ax.tick_params(length=0)
    ax.set_xlabel("recuadro rojo = élite  ·  B = banco", fontsize=7, labelpad=2)
    t = "Genes de la población"
    if generacion is not None:
        t += f", generación {generacion}"
    ax.set_title(f"{t} ({filas} mejores de {total})", fontsize=8.5)


# ---------------------------------------------------------------------------
# Matriz de rendimiento
# ---------------------------------------------------------------------------

def dibujar_matriz(ax, problema, cromosoma):
    """Rendimiento de cada jugador (filas) en cada puesto (columnas); recuadro = asignado."""
    ax.clear()
    R = problema.R
    ax.imshow(R, cmap="YlGn", aspect="auto", vmin=40, vmax=85)
    c = np.asarray(cromosoma)
    for p in range(pl.N_TITULARES):
        ax.add_patch(Rectangle((p - 0.5, c[p] - 0.5), 1, 1, fill=False, ec=ROJO, lw=1.6))
    ax.set_xticks(range(R.shape[1]))
    ax.set_xticklabels([p[0] for p in problema.puestos], fontsize=6)
    ax.set_yticks(range(R.shape[0]))
    ax.set_yticklabels([f"{j + 1} {pl.CORTOS[j]}" for j in range(R.shape[0])], fontsize=5.8)
    ax.tick_params(length=0)
    ax.set_title("Rendimiento jugador × puesto (recuadro = asignado)", fontsize=8.5)


# ---------------------------------------------------------------------------
# Ventana de la plantilla (los datos que entran a la evaluación)
# ---------------------------------------------------------------------------

def fig_plantilla():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 6.2), gridspec_kw={"width_ratios": [1.6, 1]})
    for a in (a1, a2):
        a.axis("off")
    filas = [[str(i + 1), pl.NOMBRES[i], pl.NATURAL[i], pl.PIE[i]]
             + [f"{v:.0f}" for v in pl.STATS[i]] + [f"{pl.ESTADO[i]:.0f} %"]
             for i in range(pl.N_JUGADORES)]
    t = a1.table(cellText=filas, colLabels=["#", "Jugador", "Pos.", "Pie"] + pl.ATRIBUTOS + ["Estado"],
                 loc="center", cellLoc="center")
    t.auto_set_font_size(False)
    t.set_fontsize(7.5)
    t.scale(1, 1.25)
    a1.set_title("Plantilla: estadísticas de cada jugador (1 a 99)", fontsize=10)
    filas = [[pl.ROL_NOMBRE[r]] + [f"{w:.2f}" if w else "·" for w in pesos]
             for r, pesos in pl.ROLES.items()]
    t = a2.table(cellText=filas, colLabels=["Rol"] + [a[:4] + "." for a in pl.ATRIBUTOS],
                 loc="center", cellLoc="center")
    t.auto_set_font_size(False)
    t.set_fontsize(7.5)
    t.scale(1, 1.4)
    a2.set_title("Peso de cada atributo según el rol del puesto", fontsize=10)
    fig.tight_layout()
    return fig
