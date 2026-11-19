"""
generar_resultados.py
---------------------
Corre todos los experimentos y guarda figuras y números en resultados/.

    python generar_resultados.py      (tarda unos seis minutos)
"""

import json
import os
from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import experimentos as ex
import graficas
import plantilla as pl
from algoritmo import AlgoritmoGenetico, Configuracion, n_elite
from evaluacion import ProblemaAlineacion, optimo_sin_quimica

CARPETA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultados")
os.makedirs(CARPETA, exist_ok=True)
PALETA = ["#2b6cb0", "#dd6b20", "#2f855a", "#9f7aea", "#c53030", "#718096", "#d69e2e", "#319795", "#b83280"]


def guardar(fig, nombre):
    fig.savefig(os.path.join(CARPETA, nombre), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  figura", nombre)


def figura_interfaz():
    from interfaz import Simulador
    sim = Simulador()
    for _ in range(60):
        sim.ag.paso()
    sim._dibujar()
    sim.fig.savefig(os.path.join(CARPETA, "interfaz.png"), dpi=110)
    plt.close(sim.fig)
    fig = graficas.fig_plantilla()
    guardar(fig, "plantilla.png")


def figura_evolucion():
    """La mejor alineación y los genes de la población en varias generaciones (semilla 7)."""
    p = ProblemaAlineacion("4-4-2")
    ag = AlgoritmoGenetico(p, replace(ex.BASE, semilla=7))
    ag.ejecutar()
    gens = [0, 10, 40, ag.gen_mejor]
    fig, ejes = plt.subplots(1, 4, figsize=(17, 7))
    for ax, g in zip(ejes, gens):
        graficas.dibujar_cancha(ax, p, ag.historial.mejores[g],
                                f"Generación {g} · puntaje {ag.historial.mejor[g]:.1f}")
    guardar(fig, "evolucion_alineacion.png")
    fig, ejes = plt.subplots(1, 4, figsize=(17, 4.6))
    for ax, g in zip(ejes, gens):
        graficas.dibujar_genes(ax, p, ag.historial.poblaciones[g], ag.historial.puntajes[g],
                               n_elite(ag.cfg), generacion=g)
        ax.set_title(f"Generación {g} (rojo = élite, B = banco)", fontsize=10)
    guardar(fig, "evolucion_genes.png")
    fig, ax = plt.subplots(figsize=(8, 4.2))
    graficas.dibujar_convergencia(ax, ag.historial, p.evaluar(optimo_sin_quimica(p)["cromosoma"])[0])
    guardar(fig, "convergencia_semilla7.png")
    return {"gen_mejor": ag.gen_mejor, "puntaje": ag.mejor_puntaje,
            "puntajes_gen": {g: ag.historial.mejor[g] for g in gens}}


def figura_poblacion(datos, objetivo):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 4.4))
    for i, (k, v) in enumerate(datos.items()):
        a1.plot(v["curva_media"], color=PALETA[i], lw=1.8, label=f"{k} cromosomas")
    a1.axhline(objetivo, color="k", ls=":", lw=1, label="Mejor alineación conocida")
    a1.set_ylim(objetivo - 60, objetivo + 5)
    a1.set_xlabel("Generación")
    a1.set_ylabel("Mejor puntaje (promedio de 30 corridas)")
    a1.set_title("Tamaño de la población: convergencia", fontsize=10)
    a1.grid(alpha=0.3)
    a1.legend(fontsize=8)
    claves = list(datos)
    x = np.arange(len(claves))
    a2.bar(x, [datos[k]["exito_pct"] for k in claves], color=PALETA[0], alpha=0.85)
    a2.set_xticks(x)
    a2.set_xticklabels(claves)
    a2.set_xlabel("Cromosomas en la población")
    a2.set_ylabel("Corridas que llegan a la mejor (%)")
    a2.set_ylim(0, 105)
    ax2 = a2.twinx()
    ax2.plot(x, [datos[k]["gen_llegada_media"] or np.nan for k in claves], "o-", color=PALETA[1])
    ax2.set_ylabel("Generación promedio de llegada", color=PALETA[1])
    a2.set_title("Tamaño de la población: éxito y velocidad", fontsize=10)
    a2.grid(axis="y", alpha=0.3)
    guardar(fig, "efecto_poblacion.png")


def figura_operadores(tabla):
    cruces = list(ex.CRUCES)
    muts = list(ex.MUTACIONES)
    M = np.array([[tabla[f"{c} + {m}"]["exito_pct"] for m in muts] for c in cruces])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 4.2))
    im = a1.imshow(M, cmap="Blues", vmin=0, vmax=100)
    for i in range(len(cruces)):
        for j in range(len(muts)):
            v = tabla[f"{cruces[i]} + {muts[j]}"]
            a1.text(j, i, f"{M[i, j]:.0f} %\n{v['puntaje_medio']:.1f}", ha="center", va="center",
                    color="white" if M[i, j] > 55 else "k", fontsize=9)
    a1.set_xticks(range(len(muts)))
    a1.set_xticklabels(muts)
    a1.set_yticks(range(len(cruces)))
    a1.set_yticklabels(cruces)
    a1.set_title("Éxito (%) y puntaje medio por cruce y mutación", fontsize=10)
    fig.colorbar(im, ax=a1, fraction=0.04)
    for i, k in enumerate(["PMX + Intercambio", "OX (orden) + Intercambio", "Ciclos (CX) + Intercambio",
                           "PMX + Inserción", "PMX + Inversión"]):
        a2.plot(tabla[k]["curva_media"], color=PALETA[i], lw=1.6, label=k)
    obj = max(max(v["curva_media"]) for v in tabla.values())
    a2.set_ylim(obj - 50, obj + 4)
    a2.set_xlabel("Generación")
    a2.set_ylabel("Mejor puntaje (promedio de 30)")
    a2.set_title("Convergencia de algunas combinaciones", fontsize=10)
    a2.grid(alpha=0.3)
    a2.legend(fontsize=7.5)
    guardar(fig, "comparacion_operadores.png")


def figura_parametros(barridos):
    fig, ejes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, (clave, titulo) in zip(ejes, (("elitismo_pct", "Elitismo (%)"),
                                          ("prob_mutacion", "Probabilidad de mutación"),
                                          ("seleccion", "Selección"))):
        d = barridos[clave]
        k = list(d)
        x = np.arange(len(k))
        ax.bar(x, [d[v]["exito_pct"] for v in k], color=PALETA[0], alpha=0.85)
        ax.set_ylim(0, 105)
        ax.set_xticks(x)
        ax.set_xticklabels(k)
        ax.set_xlabel(titulo)
        ax.set_ylabel("Éxito (%)")
        a2 = ax.twinx()
        a2.plot(x, [d[v]["puntaje_medio"] for v in k], "o-", color=PALETA[4])
        a2.set_ylabel("Puntaje medio", color=PALETA[4])
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Efecto de los parámetros (4-4-2, población 50, 300 generaciones, 30 semillas)", fontsize=11)
    fig.tight_layout()
    guardar(fig, "efecto_parametros.png")


def figura_sin_elitismo():
    """Sin elitismo el mejor de cada generación puede empeorar; con elitismo nunca baja."""
    p = ProblemaAlineacion("4-4-2")
    fig, ax = plt.subplots(figsize=(8, 4))
    for pct, color in ((0, PALETA[4]), (10, PALETA[0])):
        ag = AlgoritmoGenetico(p, replace(ex.BASE, elitismo_pct=pct, generaciones=150, semilla=7), False)
        ag.ejecutar()
        ax.plot(ag.historial.mejor, color=color, lw=1.4,
                label=f"Elitismo {pct} %: mejor de cada generación")
    ax.set_xlabel("Generación")
    ax.set_ylabel("Puntaje")
    ax.set_ylim(780, 860)
    ax.set_title("El elitismo impide que la mejor solución se pierda", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    guardar(fig, "efecto_elitismo.png")


def figura_formaciones(datos, lesion_datos):
    fig, ejes = plt.subplots(1, 4, figsize=(18, 7))
    for ax, f in zip(ejes, pl.FORMACIONES):
        p = ProblemaAlineacion(f)
        v, c = ex.mejor_conocido(p)
        graficas.dibujar_cancha(ax, p, c, f"{f} · puntaje {v:.1f}")
    p = ProblemaAlineacion("4-4-2", lesionados=[pl.CORTOS.index("Restrepo")])
    v, c = ex.mejor_conocido(p)
    graficas.dibujar_cancha(ejes[3], p, c, f"4-4-2 sin Restrepo (lesionado) · {v:.1f}")
    guardar(fig, "formaciones_y_lesion.png")


def figura_importancia(datos):
    """Cuánto cae el puntaje si se lesiona cada titular, y quién lo reemplaza."""
    js = datos["jugadores"][::-1]
    fig, ax = plt.subplots(figsize=(10, 5))
    y = np.arange(len(js))
    ax.barh(y, [d["caida"] for d in js], color=PALETA[4], alpha=0.85)
    for yi, d in zip(y, js):
        cambio = ", ".join(d["entran"] + d["se_mueven"])
        ax.text(d["caida"] + 0.6, yi, cambio, va="center", fontsize=7.5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{d['corto']} ({d['puesto']})" for d in js], fontsize=8.5)
    ax.set_xlabel("Puntos que pierde el equipo si ese titular no puede jugar")
    ax.set_xlim(0, max(d["caida"] for d in js) * 1.9)
    ax.set_title(f"Importancia de cada titular (4-4-2, alineación completa {datos['puntaje_completo']:.1f})",
                 fontsize=10)
    ax.grid(axis="x", alpha=0.3)
    guardar(fig, "importancia_titulares.png")


def figura_quimica(datos):
    """Rendimiento individual contra química según el peso λ."""
    lams = list(datos)
    x = np.arange(len(lams))
    rend = [datos[l]["rendimiento"] for l in lams]
    quim = [datos[l]["quimica"] for l in lams]
    fig, ax = plt.subplots(figsize=(10, 4.6))
    ax.bar(x, rend, color=PALETA[0], alpha=0.85, label="Rendimiento de los titulares")
    ax.bar(x, quim, bottom=rend, color=PALETA[2], alpha=0.85, label="Química (λ · Σ química)")
    for xi, l in zip(x, lams):
        entra = datos[l]["cambian_frente_a_lambda0"]
        texto = "Entran: " + ", ".join(entra) if entra else "Alineación sin química"
        ax.text(xi, rend[xi] + quim[xi] + 3, f"{rend[xi] + quim[xi]:.1f}\n{texto}", ha="center",
                va="bottom", fontsize=7.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"λ = {l}" for l in lams])
    ax.set_ylim(780, max(r + q for r, q in zip(rend, quim)) + 45)
    ax.set_ylabel("Puntaje del equipo")
    ax.set_title("Peso de la química: cuándo cambia la mejor alineación (4-4-2)", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(axis="y", alpha=0.3)
    guardar(fig, "efecto_quimica.png")


def convertir(o):
    if isinstance(o, dict):
        return {str(k): convertir(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [convertir(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def principal():
    p = ProblemaAlineacion("4-4-2")
    objetivo, _ = ex.mejor_conocido(p)
    print("Interfaz y plantilla...")
    figura_interfaz()
    print("Evolución de una corrida...")
    evol = figura_evolucion()
    print("Tamaño de la población...")
    pob = ex.barrido(p, "tam_poblacion", [10, 20, 50, 100, 200])
    figura_poblacion(pob, objetivo)
    print("Operadores...")
    ops = ex.operadores(p)
    figura_operadores(ops)
    print("Elitismo, mutación y selección...")
    barridos = {"elitismo_pct": ex.barrido(p, "elitismo_pct", [0, 2, 10, 20, 40]),
                "prob_mutacion": ex.barrido(p, "prob_mutacion", [0.0, 0.05, 0.1, 0.2, 0.4]),
                "seleccion": ex.selecciones(p)}
    figura_parametros(barridos)
    figura_sin_elitismo()
    print("Química contra método húngaro, formaciones y lesión...")
    quim = ex.quimica_contra_hungaro()
    form = ex.formaciones()
    les = ex.lesion()
    figura_formaciones(form, les)
    print("Experimentos con el fenómeno: importancia de cada titular y peso de la química...")
    imp = ex.importancia_jugadores()
    figura_importancia(imp)
    bq = ex.barrido_quimica()
    figura_quimica(bq)

    for grupo in [pob, ops] + list(barridos.values()):
        for v in grupo.values():
            v.pop("curva_media", None)
            v.pop("curva_media_poblacion", None)
    with open(os.path.join(CARPETA, "resultados.json"), "w", encoding="utf-8") as fh:
        json.dump(convertir({"mejor_conocido_442": objetivo, "evolucion_semilla7": evol,
                             "poblacion": pob, "operadores": ops, "parametros": barridos,
                             "quimica_vs_hungaro": quim, "formaciones": form, "lesion": les,
                             "importancia_titulares": imp, "barrido_quimica": bq}),
                  fh, ensure_ascii=False, indent=1)
    print("Listo: resultados/resultados.json")


if __name__ == "__main__":
    principal()
