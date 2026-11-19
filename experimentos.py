"""
experimentos.py
---------------
Experimentos repetidos con 30 semillas: como el algoritmo es aleatorio,
cada corrida puede dar algo distinto, así que se reporta el promedio, la
desviación y el porcentaje de corridas que llegan a la mejor alineación
conocida.
"""

from dataclasses import replace

import numpy as np

import plantilla as pl
from algoritmo import AlgoritmoGenetico, Configuracion
from evaluacion import ProblemaAlineacion, optimo_sin_quimica
from operadores import CRUCES, MUTACIONES, SELECCIONES

SEMILLAS = list(range(1, 31))
BASE = Configuracion(tam_poblacion=50, generaciones=300)
_mejor_conocido = {}


def mejor_conocido(problema, corridas=8):
    """Mejor alineación encontrada con corridas largas (población 300, 400 generaciones).

    Con química no hay un método exacto que la garantice, así que se usa como
    referencia la mejor de varias corridas largas y se comprueba que ningún
    intercambio de dos jugadores la mejora.
    """
    clave = (problema.formacion, problema.lambda_quimica, tuple(problema.lesionados.tolist()))
    if clave not in _mejor_conocido:
        mejor = None
        for s in range(corridas):
            r = AlgoritmoGenetico(problema, Configuracion(tam_poblacion=300, generaciones=400,
                                                          semilla=1000 + s), False).ejecutar()
            if mejor is None or r["mejor_puntaje"] > mejor["mejor_puntaje"]:
                mejor = r
        _mejor_conocido[clave] = (mejor["mejor_puntaje"], mejor["mejor"])
    return _mejor_conocido[clave]


def es_optimo_local(problema, cromosoma):
    """True si ningún intercambio de dos genes mejora el puntaje."""
    c = np.asarray(cromosoma)
    v = problema.evaluar(c)[0]
    for i in range(c.size):
        for j in range(i + 1, c.size):
            d = c.copy()
            d[i], d[j] = d[j], d[i]
            if problema.evaluar(d)[0] > v + 1e-9:
                return False
    return True


def correr(problema, cfg, semillas=SEMILLAS):
    salida = []
    for s in semillas:
        ag = AlgoritmoGenetico(problema, replace(cfg, semilla=s), guardar_poblaciones=False)
        r = ag.ejecutar()
        salida.append((r, ag.historial))
    return salida


def resumen(problema, cfg, semillas=SEMILLAS, objetivo=None):
    """Estadísticas de una configuración sobre 30 semillas."""
    if objetivo is None:
        objetivo = mejor_conocido(problema)[0]
    corridas = correr(problema, cfg, semillas)
    finales = np.array([r["mejor_puntaje"] for r, _ in corridas])
    exito = finales >= objetivo - 1e-6
    gen_llegada = [int(np.argmax(np.array(h.mejor_global) >= objetivo - 1e-6))
                   for (r, h), ok in zip(corridas, exito) if ok]
    curvas = np.array([h.mejor_global for _, h in corridas])
    distintas = len({tuple(r["mejor"][:pl.N_TITULARES]) for r, _ in corridas})
    return {
        "puntaje_medio": float(finales.mean()), "puntaje_desv": float(finales.std()),
        "puntaje_peor": float(finales.min()), "puntaje_mejor": float(finales.max()),
        "exito_pct": float(100 * exito.mean()),
        "gen_llegada_media": float(np.mean(gen_llegada)) if gen_llegada else None,
        "gen_llegada_max": int(np.max(gen_llegada)) if gen_llegada else None,
        "evaluaciones_llegada": float(np.mean(gen_llegada) * cfg.tam_poblacion) if gen_llegada else None,
        "alineaciones_distintas": distintas,
        "diversidad_final": float(np.mean([r["diversidad_final"] for r, _ in corridas])),
        "curva_media": curvas.mean(axis=0).tolist(),
        "curva_media_poblacion": np.mean([h.media for _, h in corridas], axis=0).tolist(),
    }


def barrido(problema, parametro, valores, base=BASE, semillas=SEMILLAS):
    return {str(v): resumen(problema, replace(base, **{parametro: v}), semillas) for v in valores}


def operadores(problema, base=BASE, semillas=SEMILLAS):
    tabla = {}
    for c in CRUCES:
        for m in MUTACIONES:
            tabla[f"{c} + {m}"] = resumen(problema, replace(base, cruce=c, mutacion=m), semillas)
    return tabla


def selecciones(problema, base=BASE, semillas=SEMILLAS):
    return {s: resumen(problema, replace(base, seleccion=s), semillas) for s in SELECCIONES}


def quimica_contra_hungaro(formacion="4-4-2", semillas=SEMILLAS):
    """Sin química el húngaro da el óptimo exacto: ¿el AG lo alcanza? Con química, ¿lo supera?"""
    salida = {}
    for lam in (0.0, 3.0):
        p = ProblemaAlineacion(formacion, lam)
        h = optimo_sin_quimica(p)
        r = resumen(p, BASE, semillas)
        mejor, crom = mejor_conocido(p)
        salida[str(lam)] = {"hungaro": float(p.evaluar(h["cromosoma"])[0]),
                            "hungaro_alineacion": [pl.CORTOS[j] for j in h["cromosoma"][:11]],
                            "ag_mejor": mejor, "ag_alineacion": [pl.CORTOS[j] for j in crom[:11]],
                            "ag_medio": r["puntaje_medio"], "exito_pct": r["exito_pct"]}
    return salida


def formaciones(lam=3.0):
    salida = {}
    for f in pl.FORMACIONES:
        p = ProblemaAlineacion(f, lam)
        v, c = mejor_conocido(p)
        d = p.detalle(c)
        salida[f] = {"puntaje": v, "rendimiento": d["rendimiento"], "quimica": d["quimica"],
                     "alineacion": [(x["puesto"], x["corto"], round(x["rendimiento"], 1), x["natural"])
                                    for x in d["titulares"]],
                     "banco": d["banco"], "optimo_local": es_optimo_local(p, c)}
    return salida


def lesion(formacion="4-4-2", lesionado="Restrepo"):
    """Se lesiona un titular: ¿cómo reacomoda el algoritmo al equipo?"""
    j = pl.CORTOS.index(lesionado)
    antes = ProblemaAlineacion(formacion)
    despues = ProblemaAlineacion(formacion, lesionados=[j])
    va, ca = mejor_conocido(antes)
    vd, cd = mejor_conocido(despues)
    # Se compara por jugador y sigla del puesto, así dos centrales que solo
    # cambian de lado no cuentan como cambio.
    pos_a = {pl.CORTOS[ca[p]]: antes.puestos[p][0] for p in range(pl.N_TITULARES)}
    pos_d = {pl.CORTOS[cd[p]]: despues.puestos[p][0] for p in range(pl.N_TITULARES)}
    cambios = [(nombre, pos_a.get(nombre, "banco"), pos_d.get(nombre, "banco"))
               for nombre in pl.CORTOS if pos_a.get(nombre) != pos_d.get(nombre)]
    return {"lesionado": pl.NOMBRES[j], "puntaje_antes": va, "puntaje_despues": vd,
            "alineacion_despues": [(despues.puestos[p][0], pl.CORTOS[cd[p]]) for p in range(11)],
            "cambios": cambios}
