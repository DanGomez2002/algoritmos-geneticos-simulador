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
    elif len(g) == 1:
        ax.set_xlim(0, 1)            # generación 0: sin esto el eje sale de -0.05 a 0.05
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


def dibujar_genes(ax, problema, poblacion, puntajes, n_elite=0, max_filas=15, generacion=None, letra=6):
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
    if letra:
        _numeros_en_celdas(ax, P, fontsize=letra)
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
    ax.set_title(f"{t} ({filas} mejores de {total})" if filas < total else f"{t} (los {total} cromosomas)",
                 fontsize=8.5)


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
    """Ventana «Ver plantilla»: jugadores arriba y pesos de cada rol debajo (cabe a lo ancho en un portátil)."""
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 9.4), gridspec_kw={"height_ratios": [1.9, 1]})
    for a in (a1, a2):
        a.axis("off")
    filas = [[str(i + 1), pl.NOMBRES[i], pl.NATURAL[i], pl.PIE[i]]
             + [f"{v:.0f}" for v in pl.STATS[i]] + [f"{pl.ESTADO[i]:.0f} %"]
             for i in range(pl.N_JUGADORES)]
    t = a1.table(cellText=filas, colLabels=["#", "Jugador", "Pos.", "Pie"] + pl.ATRIBUTOS + ["Estado"],
                 loc="center", cellLoc="center")
    t.auto_set_font_size(False)
    t.set_fontsize(8)
    t.scale(1, 1.2)
    t.auto_set_column_width([1])            # que el nombre completo quepa en su celda
    a1.set_title("Plantilla: estadísticas de cada jugador (1 a 99)", fontsize=10)
    filas = [[pl.ROL_NOMBRE[r]] + [f"{w:.2f}" if w else "·" for w in pesos]
             for r, pesos in pl.ROLES.items()]
    t = a2.table(cellText=filas, colLabels=["Rol"] + pl.ATRIBUTOS, loc="center", cellLoc="center")
    t.auto_set_font_size(False)
    t.set_fontsize(8)
    t.scale(1, 1.2)
    t.auto_set_column_width([0])
    a2.set_title("Peso de cada atributo según el rol del puesto (cada fila suma 1)", fontsize=10)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Ventana «Ver genes»: la tabla de genes en grande
# ---------------------------------------------------------------------------

def generaciones_con_mejora(historial):
    """Generaciones en que la mejor alineación encontrada mejoró (la 0 incluida)."""
    g = historial.mejor_global
    return [0] + [k for k in range(1, len(g)) if g[k] > g[k - 1] + 1e-9]


def letra_para_filas(ax, filas, maxima=7.0, minima=3.2):
    """Tamaño de los números que cabe en filas de la altura actual del eje (0 si no se leerían)."""
    alto_pt = ax.get_position().height * ax.figure.get_figheight() * 72 / max(1, filas)
    letra = min(maxima, 0.7 * alto_pt)
    return letra if letra >= minima else 0


def dibujar_mejores(ax, problema, historial, hasta, max_filas=30):
    """Una fila por cada generación en que mejoró la mejor alineación, hasta la generación `hasta`.

    Muestra cómo se fueron acomodando los genes del mejor cromosoma: las
    primeras filas son mosaicos distintos y las últimas solo cambian uno o
    dos jugadores.
    """
    ax.clear()
    gens = [k for k in generaciones_con_mejora(historial) if k <= hasta]
    if len(gens) > max_filas:                       # la primera y las más recientes
        gens = gens[:1] + gens[-(max_filas - 1):]
    P = np.array([historial.mejores[k] for k in gens])
    imagen = np.array([[COLORES_JUGADOR[j] for j in fila] for fila in P])
    ax.imshow(imagen, aspect="auto", interpolation="nearest")
    letra = letra_para_filas(ax, len(gens), maxima=8.5)
    if letra:
        _numeros_en_celdas(ax, P, fontsize=letra)
    ax.axvline(pl.N_TITULARES - 0.5, color="k", lw=1.5)
    L = P.shape[1]
    siglas = [p[0] for p in problema.puestos] + ["B"] * (L - pl.N_TITULARES)
    ax.set_xticks(range(L))
    ax.set_xticklabels(siglas, fontsize=8)
    ax.set_yticks(range(len(gens)))
    ax.set_yticklabels([f"gen. {k} · {historial.mejor_global[k]:.1f}" for k in gens], fontsize=8)
    ax.tick_params(length=0)
    ax.set_xlabel("B = banco  ·  cada número es un jugador (el de «Ver plantilla»)", fontsize=7.5, labelpad=3)
    ax.set_title("Mejor cromosoma cada vez que mejoró la alineación (generación · puntaje)", fontsize=10)


# ---------------------------------------------------------------------------
# Ventana «Resumen»: cómo avanzó la búsqueda hasta la generación actual
# ---------------------------------------------------------------------------

def cambios_entre(problema, antes, despues):
    """Qué cambió en los titulares entre dos alineaciones, en pocas palabras.

    «+Murcia LI» entra desde el banco, «Silva MC→DFC» cambia de puesto y
    «−Calderón» sale al banco. Pasar de un puesto a otro igual (los dos
    centrales o los dos delanteros que cambian de lado) no se cuenta.
    """
    n = pl.N_TITULARES
    siglas = [p[0] for p in problema.puestos]
    puesto_antes = {int(j): k for k, j in enumerate(antes[:n])}
    partes = []
    for k in range(n):
        j = int(despues[k])
        if int(antes[k]) == j:
            continue
        if j in puesto_antes:
            if siglas[puesto_antes[j]] == siglas[k]:
                continue
            partes.append(f"{pl.CORTOS[j]} {siglas[puesto_antes[j]]}→{siglas[k]}")
        else:
            partes.append(f"+{pl.CORTOS[j]} {siglas[k]}")
    titulares_despues = set(int(j) for j in despues[:n])
    partes += [f"−{pl.CORTOS[j]}" for j in antes[:n] if int(j) not in titulares_despues]
    return " · ".join(partes) if partes else "solo cambia el banco"


def fig_resumen(ag, referencia=None, max_filas=16):
    """Resumen de la corrida: evolución, diversidad, mejoras generación a generación y datos finales."""
    h, cfg, problema = ag.historial, ag.cfg, ag.problema
    fig = plt.figure(figsize=(14, 8.6))
    fig.suptitle(f"Resumen de la corrida · formación {problema.formacion} · generación {ag.generacion}",
                 fontsize=12, weight="bold", x=0.01, ha="left")
    ax_ev = fig.add_axes([0.05, 0.56, 0.56, 0.36])
    ax_div = fig.add_axes([0.68, 0.56, 0.30, 0.36])
    ax_tab = fig.add_axes([0.01, 0.07, 0.62, 0.38])
    ax_txt = fig.add_axes([0.66, 0.07, 0.33, 0.41])

    gens = np.arange(len(h.mejor))
    ax_ev.fill_between(gens, h.peor, h.mejor, color=GRIS, alpha=0.15, lw=0, label="Peor a mejor de la generación")
    ax_ev.plot(gens, h.media, color=GRIS, lw=1.2, label="Promedio de la población")
    ax_ev.plot(gens, h.mejor_global, color=ROJO, lw=2.2, label="Mejor encontrado (élite)")
    mejoras = generaciones_con_mejora(h)
    ax_ev.plot(mejoras, [h.mejor_global[k] for k in mejoras], "o", color=ROJO, ms=4)
    if referencia is not None:
        ax_ev.axhline(referencia, color=VERDE, ls="--", lw=1.2, label="Asignación húngara (ignora la química)")
    bajo = min(np.percentile(h.media, 5), referencia if referencia is not None else np.inf) - 15
    ax_ev.set_ylim(bajo, max(h.mejor_global) + 8)
    ax_ev.set_xlim(0, max(1, len(gens) - 1))
    ax_ev.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax_ev.set_title("Evolución del puntaje (los puntos marcan cada mejora)", fontsize=10)
    ax_ev.set_xlabel("Generación", fontsize=9)
    ax_ev.set_ylabel("Puntaje del equipo", fontsize=9)
    ax_ev.grid(alpha=0.3)
    ax_ev.legend(fontsize=7.5, loc="lower right")

    ax_div.plot(gens, h.diversidad, color=AZUL, lw=1.5)
    ax_div.set_ylim(0, 100)
    ax_div.set_xlim(0, max(1, len(gens) - 1))
    ax_div.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax_div.set_title("Diversidad: genes distintos del mejor (%)", fontsize=10)
    ax_div.set_xlabel("Generación", fontsize=9)
    ax_div.grid(alpha=0.3)

    ax_tab.axis("off")
    filas = []
    for i, g in enumerate(mejoras):
        if i == 0:
            que, delta = "población inicial (permutaciones al azar)", ""
        else:
            que = cambios_entre(problema, h.mejores[mejoras[i - 1]], h.mejores[g])
            delta = f"+{h.mejor_global[g] - h.mejor_global[mejoras[i - 1]]:.1f}"
        if len(que) > 72:
            que = que[:71] + "…"
        filas.append([str(g), f"{h.mejor_global[g]:.1f}", delta, que])
    if len(filas) > max_filas:                      # la primera y las últimas
        omitidas = len(filas) - max_filas + 1
        filas = filas[:1] + [["…", "", "", f"{omitidas} mejoras más"]] + filas[-(max_filas - 2):]
    t = ax_tab.table(cellText=filas, colLabels=["Gen.", "Puntaje", "Mejora", "Qué cambió en los titulares"],
                     colWidths=[0.07, 0.1, 0.09, 0.74], loc="upper center", cellLoc="center")
    t.auto_set_font_size(False)
    t.set_fontsize(7.8)
    t.scale(1, 1.18)
    for (r, c), celda in t.get_celld().items():
        celda.set_edgecolor("#cbd5e0")
        if r == 0:
            celda.set_facecolor("#edf2f7")
            celda.set_text_props(weight="bold")
        elif c == 3:
            celda.set_text_props(ha="left")
            celda.PAD = 0.01
    ax_tab.set_title("Generaciones en que mejoró la mejor alineación", fontsize=10)

    ax_txt.axis("off")
    d = problema.detalle(ag.mejor)
    parada = [f"máx. {cfg.generaciones}"]
    if cfg.sin_mejora:
        parada.append(f"{cfg.sin_mejora} sin mejora")
    if cfg.umbral > 0:
        parada.append(f"objetivo {cfg.umbral:g}")
    lineas = [
        "CONFIGURACIÓN",
        f"Población {cfg.tam_poblacion} · semilla {cfg.semilla} · λ = {problema.lambda_quimica:g}",
        f"{cfg.seleccion} · {cfg.cruce} · {cfg.mutacion}",
        f"pc {cfg.prob_cruce:.2f} · pm {cfg.prob_mutacion:.2f} · élite {cfg.elitismo_pct:g} %",
        "Parada: " + " · ".join(parada),
        "",
        "RESULTADO",
        f"Generaciones corridas  {ag.generacion}",
        f"Motivo   {ag.motivo_parada() or 'la corrida sigue'}",
        f"Mejor puntaje  {d['total']:.2f} (gen. {ag.gen_mejor})",
        f"  rendimiento {d['rendimiento']:.2f} · química {d['quimica']:+.0f}",
        f"  penalización {-d['penalizacion'] + 0.0:+.0f}",
    ]
    if referencia is not None:
        lineas.append(f"Frente al húngaro  {d['total'] - referencia:+.2f}")
    lineas += [
        f"Mejoras  {len(mejoras) - 1} · evaluaciones {ag.evaluaciones:,}".replace(",", "."),
        f"Diversidad final  {h.diversidad[-1]:.1f} %",
        "",
        "MEJOR ALINEACIÓN",
    ]
    tit = d["titulares"]
    for k in range(0, len(tit), 2):
        lineas.append("   ".join(f"{f['puesto']:<4}{f['corto']:<9}" for f in tit[k:k + 2]))
    lineas += textwrap.wrap("Banco: " + ", ".join(pl.CORTOS[pl.NOMBRES.index(b)] for b in d["banco"]), 46,
                            subsequent_indent="       ")
    ax_txt.text(0, 1, "\n".join(lineas), va="top", ha="left", fontsize=8.2, family="monospace",
                transform=ax_txt.transAxes)
    return fig
