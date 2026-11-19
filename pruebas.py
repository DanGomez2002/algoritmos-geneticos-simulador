"""
pruebas.py
----------
Comprobaciones automáticas. Correr con:  python pruebas.py
"""

import os
import tempfile
import traceback
from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import numpy as np

import abrir_codigo
import plantilla as pl
from algoritmo import AlgoritmoGenetico, Configuracion, aptitud_desde_puntaje, n_elite
from evaluacion import ProblemaAlineacion, matriz_rendimiento, optimo_sin_quimica, rendimiento_base
from operadores import (CRUCES, MUTACIONES, ciclos, es_permutacion, ox, pmx,
                        seleccion_ruleta, seleccion_torneo, seleccion_ranking)

RNG = lambda s=0: np.random.default_rng(s)


def test_datos_consistentes():
    assert pl.STATS.shape == (18, 7)
    for rol, pesos in pl.ROLES.items():
        assert abs(sum(pesos) - 1) < 1e-9, rol
    for f, puestos in pl.FORMACIONES.items():
        assert len(puestos) == 11 and puestos[0][1] == "POR", f
    assert np.allclose(pl.QUIMICA, pl.QUIMICA.T)


def test_rendimiento_a_mano():
    # Ríos (portero): 0,05·45 + 0,05·55 + 0,05·70 + 0,85·82 = 78,2
    assert abs(rendimiento_base(0, "POR") - 78.2) < 1e-9
    R = matriz_rendimiento("4-4-2")
    # Murcia es zurdo: de lateral derecho pierde 8 % frente a la misma suma ponderada.
    base = rendimiento_base(6, "LAT") * (0.7 + 0.3 * 0.90)
    assert abs(R[6, 4] - base) < 1e-9 and abs(R[6, 1] - 0.92 * base) < 1e-9


def test_cruces_dan_permutaciones():
    rng = RNG()
    for _ in range(500):
        a, b = rng.permutation(18), rng.permutation(18)
        for nombre, f in CRUCES.items():
            h1, h2 = f(a, b, rng)
            assert es_permutacion(h1, 18) and es_permutacion(h2, 18), nombre
        for nombre, m in MUTACIONES.items():
            assert es_permutacion(m(a, rng), 18), nombre


def test_pmx_ejemplo_de_libro():
    # Cortes en 3 y 7 (posiciones 3 a 6), resuelto a mano.
    p1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]) - 1
    p2 = np.array([4, 5, 2, 1, 8, 7, 6, 9, 3]) - 1
    h1, _ = pmx(p1, p2, None, cortes=(3, 7))
    # Tramo copiado de p1: 4 5 6 7; mapa 4->1, 5->8, 6->7, 7->6.
    assert (h1 + 1).tolist() == [1, 8, 2, 4, 5, 6, 7, 9, 3]


def test_ox_ejemplo_de_libro():
    p1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9]) - 1
    p2 = np.array([9, 3, 7, 8, 2, 6, 5, 1, 4]) - 1
    h1, _ = ox(p1, p2, None, cortes=(3, 7))
    assert (h1 + 1).tolist() == [3, 8, 2, 4, 5, 6, 7, 1, 9]


def test_ciclos_ejemplo_de_libro():
    p1 = np.array([1, 2, 3, 4, 5, 6, 7, 8]) - 1
    p2 = np.array([8, 5, 2, 1, 3, 6, 4, 7]) - 1
    h1, h2 = ciclos(p1, p2)
    assert (h1 + 1).tolist() == [1, 5, 2, 4, 3, 6, 7, 8]
    # Cada gen del hijo está en una posición que tenía en alguno de los padres.
    assert np.all((h1 == p1) | (h1 == p2))


def test_selecciones():
    apt = np.array([1.0, 3.0, 6.0])
    frec = np.bincount(seleccion_ruleta(apt, 60000, RNG()), minlength=3) / 60000
    assert np.allclose(frec, [0.1, 0.3, 0.6], atol=0.01)
    for sel in (seleccion_torneo, seleccion_ranking):
        assert np.mean(sel(np.arange(10.0), 20000, RNG())) > 5


def test_aptitud_todos_pueden_ser_padres():
    a = aptitud_desde_puntaje([800.0, 700.0, 750.0])
    assert a.argmax() == 0 and a.min() > 0


def test_elitismo():
    assert n_elite(Configuracion(tam_poblacion=50, elitismo_pct=10)) == 5
    assert n_elite(Configuracion(tam_poblacion=10, elitismo_pct=10)) == 1
    assert n_elite(Configuracion(tam_poblacion=10, elitismo_pct=1)) == 1
    assert n_elite(Configuracion(elitismo_pct=0)) == 0
    ag = AlgoritmoGenetico(ProblemaAlineacion(), Configuracion(generaciones=80, semilla=3))
    ag.ejecutar()
    assert np.all(np.diff(ag.historial.mejor) >= -1e-9)        # con élite nunca baja


def test_poblacion_siempre_valida():
    ag = AlgoritmoGenetico(ProblemaAlineacion("3-5-2"), Configuracion(generaciones=30, cruce="OX (orden)",
                                                                      mutacion="Inversión"))
    ag.ejecutar()
    for P in ag.historial.poblaciones:
        assert all(es_permutacion(c, 18) for c in P)


def test_penalizaciones():
    p = ProblemaAlineacion(lesionados=[9])
    c = np.arange(18)                       # Restrepo (9) titular; porteros 0 y 1 en cancha
    d = p.detalle(c)
    assert d["penalizacion"] == 300 + 40
    c2 = c.copy()
    c2[[1, 11]] = c2[[11, 1]]               # Cuéllar al banco
    c2[[9, 12]] = c2[[12, 9]]               # Restrepo al banco
    assert p.detalle(c2)["penalizacion"] == 0


def test_portero_lesionado_no_es_suplente():
    p = ProblemaAlineacion(lesionados=[1])
    c = np.arange(18)
    c[[1, 11]] = c[[11, 1]]                 # Cuéllar (lesionado) en el banco, Ríos de titular
    assert p.detalle(c)["penalizacion"] == 40
    sano = ProblemaAlineacion()
    assert sano.detalle(c)["penalizacion"] == 0


def test_inversion_siempre_cambia():
    rng = RNG()
    for _ in range(2000):
        c = rng.permutation(18)
        assert not np.array_equal(MUTACIONES["Inversión"](c, rng), c)


def test_portero_lesionado_entra_el_suplente():
    p = ProblemaAlineacion(lesionados=[pl.CORTOS.index("Ríos")])
    r = AlgoritmoGenetico(p, Configuracion(tam_poblacion=100, generaciones=150), False).ejecutar()
    assert r["mejor"][0] == pl.CORTOS.index("Cuéllar")


def test_extremo_de_volante_no_esta_fuera_de_posicion():
    import graficas
    monroy = pl.CORTOS.index("Monroy")
    assert not graficas.fuera_de_posicion(monroy, "MED") and not graficas.fuera_de_posicion(monroy, "EXT")
    assert graficas.fuera_de_posicion(monroy, "DFC")


def test_hungaro_sin_quimica_lo_alcanza_el_ag():
    p = ProblemaAlineacion(lambda_quimica=0)
    h = optimo_sin_quimica(p)
    r = AlgoritmoGenetico(p, Configuracion(tam_poblacion=200, generaciones=300, semilla=1), False).ejecutar()
    assert abs(r["mejor_puntaje"] - h["rendimiento"]) < 1e-6


def test_con_quimica_supera_al_hungaro():
    p = ProblemaAlineacion(lambda_quimica=3)
    h = p.evaluar(optimo_sin_quimica(p)["cromosoma"])[0]
    r = AlgoritmoGenetico(p, Configuracion(tam_poblacion=100, generaciones=300, semilla=2), False).ejecutar()
    assert r["mejor_puntaje"] > h + 5


def test_lesionado_no_juega():
    j = pl.CORTOS.index("Restrepo")
    p = ProblemaAlineacion(lesionados=[j])
    r = AlgoritmoGenetico(p, Configuracion(generaciones=150), False).ejecutar()
    assert j not in r["mejor"][:11]


def test_reproducible_y_parada():
    a = AlgoritmoGenetico(ProblemaAlineacion(), Configuracion(semilla=5)).ejecutar()
    b = AlgoritmoGenetico(ProblemaAlineacion(), Configuracion(semilla=5)).ejecutar()
    assert a["mejor_puntaje"] == b["mejor_puntaje"]
    ag = AlgoritmoGenetico(ProblemaAlineacion(), Configuracion(generaciones=1000, sin_mejora=20))
    ag.ejecutar()
    assert ag.generacion < 1000 and "sin mejorar" in ag.motivo_parada()
    ag = AlgoritmoGenetico(ProblemaAlineacion(), Configuracion(generaciones=1000, umbral=800))
    ag.ejecutar()
    assert ag.mejor_puntaje >= 800 and "objetivo" in ag.motivo_parada()


# ------------------------------------------------------------------ interfaz
def _evento(sim, ax, boton=3, y=None, key=None):
    from matplotlib.backend_bases import MouseEvent
    sim.fig.canvas.draw()
    caja = ax.get_window_extent()
    ev = MouseEvent("button_press_event", sim.fig.canvas, caja.x0 + caja.width / 2,
                    caja.y0 + caja.height / 2 if y is None else y, button=boton, key=key)
    ev.inaxes = ax
    return ev


def test_interfaz_controles():
    from interfaz import Simulador
    sim = Simulador()
    for _ in range(4):
        sim._un_paso()
    assert sim.ag.generacion == 4
    sim.s_pm.set_val(0.3)                                    # en vivo: no reinicia
    assert sim.ag.generacion == 4 and sim.ag.cfg.prob_mutacion == 0.3
    sim.r_cru.set_active(2)
    assert sim.ag.cfg.cruce == "Ciclos (CX)" and sim.ag.generacion == 4
    sim._ver_generacion(1)                                   # volver a una generación vieja
    assert sim.gen_vista == 1
    sim._elegir_formacion("4-3-3")
    assert sim.ag.generacion == 0 and sim.problema.formacion == "4-3-3"
    sim.t_les.set_val("10, 3")
    assert set(sim.problema.lesionados.tolist()) == {9, 2}
    sim.t_sem.set_val("-3")                                  # semilla negativa: no debe fallar
    assert sim.ag.cfg.semilla == 3
    sim._repetir()
    for _ in range(200):
        sim._tic()
    assert sim.gen_vista is None and not sim.repitiendo
    ruta = os.path.join(tempfile.gettempdir(), "captura_futbol.png")
    sim.fig.savefig(ruta)
    assert os.path.getsize(ruta) > 10000


def test_clic_derecho():
    from interfaz import Simulador
    pedidos = []
    original = abrir_codigo.abrir
    abrir_codigo.abrir = lambda a, n, b=None: (pedidos.append((a, n)) or (True, f"{a}:1"))
    try:
        sim = Simulador()
        for ax, (archivo, nombre, otros) in sim.destinos.items():
            for a, n in [(archivo, nombre)] + list(otros):
                assert abrir_codigo.linea_de(a, n) > 1, (a, n)
        sim._clic_derecho(_evento(sim, sim.ax_can))
        assert pedidos[-1] == ("graficas.py", "dibujar_cancha")
        etiqueta = [t for t in sim.r_cru.labels if t.get_text() == "OX (orden)"][0]
        caja = etiqueta.get_window_extent()
        sim._clic_derecho(_evento(sim, sim.ax_cru, y=(caja.y0 + caja.y1) / 2))
        assert pedidos[-1] == ("operadores.py", "ox")
        sim._clic_derecho(_evento(sim, sim.ax_eli, key="shift"))
        assert pedidos[-1] == ("algoritmo.py", "AlgoritmoGenetico.paso")
        n = len(pedidos)
        sim._clic_derecho(_evento(sim, sim.ax_k, boton=1))
        assert len(pedidos) == n
        # El clic derecho sobre "Un paso" no avanza la generación.
        from matplotlib.backend_bases import MouseEvent
        antes = sim.ag.generacion
        for tipo in ("button_press_event", "button_release_event"):
            ev = _evento(sim, sim.ax_paso)
            ev.name = tipo
            sim.fig.canvas.callbacks.process(tipo, ev)
        assert sim.ag.generacion == antes
    finally:
        abrir_codigo.abrir = original


if __name__ == "__main__":
    pruebas = [(n, f) for n, f in list(globals().items()) if n.startswith("test_")]
    fallas = 0
    for nombre, f in pruebas:
        try:
            f()
            print(f"OK     {nombre}")
        except Exception:
            fallas += 1
            print(f"FALLA  {nombre}")
            traceback.print_exc()
    print(f"\n{len(pruebas) - fallas} de {len(pruebas)} comprobaciones pasan.")
    raise SystemExit(1 if fallas else 0)
