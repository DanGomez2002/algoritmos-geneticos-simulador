"""
interfaz.py
-----------
Ventana del simulador (widgets de Matplotlib).

Izquierda : formación y operadores (selección, cruce, mutación).
Centro    : parámetros, condiciones de parada, escenario (lesionados),
            botones y panel de resultados.
Derecha   : la cancha con la mejor alineación, la evolución del puntaje,
            los genes de la población y la matriz de rendimiento.

Debajo de la cancha, el deslizador "Ver generación" permite volver a
cualquier generación y ver cómo estaban organizados los genes; el botón
"Repetir evolución" la reproduce desde la generación 0.

Clic derecho sobre cualquier control o gráfica: abre su código en VS Code.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, RadioButtons, Slider, TextBox

import abrir_codigo
import graficas
import plantilla as pl
from algoritmo import AlgoritmoGenetico, Configuracion, n_elite
from evaluacion import ProblemaAlineacion, optimo_sin_quimica
from operadores import CRUCES, MUTACIONES, SELECCIONES

# Valores por defecto de la ventana (las corridas de experimentos.py usan los suyos).
POBLACION_VENTANA = 100
SIN_MEJORA_VENTANA = 50
LAMBDA_VENTANA = 3.0
# Mejor puntaje conocido de cada formación con λ = 3 y sin lesionados
# (resultados/resultados.json), redondeado hacia abajo: es el puntaje objetivo
# por defecto. Con otro λ o con lesionados el objetivo se desconoce y se apaga (0).
OBJETIVO_CONOCIDO = {"4-4-2": 851.2, "4-3-3": 848.2, "3-5-2": 845.5}


class Simulador:

    def __init__(self):
        self.formacion = "4-4-2"
        self.corriendo = False
        self.repitiendo = False
        self.gen_vista = None               # None = seguir la última generación

        self.fig = plt.figure(figsize=(16, 9))
        try:
            self.fig.canvas.manager.set_window_title("Algoritmo genético: alineación de un equipo de fútbol")
        except Exception:
            pass
        self.fig.text(0.01, 0.977, "Alineación de un equipo de fútbol con algoritmo genético",
                      fontsize=13.5, weight="bold", va="top")
        self.fig.text(0.01, 0.948, "Computacional I · Carlos Daniel Gómez Murcia", fontsize=8.5, color="#4a5568", va="top")
        self.fig.text(0.99, 0.006, "Clic derecho sobre cualquier control o gráfica: abre su código en VS Code",
                      fontsize=7.5, color="#a0aec0", ha="right")
        self.estado = self.fig.text(0.15, 0.006, "", fontsize=8, color=graficas.AZUL)

        self._crear_controles()
        self._crear_graficas()
        self.reiniciar()

        self.timer = self.fig.canvas.new_timer(interval=60)
        self.timer.add_callback(self._tic)
        self.fig.canvas.mpl_connect("button_press_event", self._clic_derecho)
        self.fig.canvas.mpl_connect("resize_event", lambda _: self._ajustar_letras())

    # ------------------------------------------------------------------ armado
    def _radio(self, rect, titulo, opciones, activa, accion):
        ax = self.fig.add_axes(rect)
        ax.set_title(titulo, fontsize=9, loc="left", weight="bold")
        r = RadioButtons(ax, opciones, active=opciones.index(activa))
        for t in r.labels:
            t.set_fontsize(8.5)
        r.on_clicked(accion)
        return ax, r

    def _slider(self, y, etiqueta, vmin, vmax, valor, paso, fmt="%d"):
        ax = self.fig.add_axes([0.215, y, 0.10, 0.02])
        s = Slider(ax, etiqueta, vmin, vmax, valinit=valor, valstep=paso, valfmt=fmt)
        s.label.set_fontsize(8)
        s.valtext.set_fontsize(8)
        return ax, s

    def _caja(self, y, etiqueta, inicial, accion):
        ax = self.fig.add_axes([0.215, y, 0.07, 0.026])
        t = TextBox(ax, etiqueta + " ", initial=inicial)
        t.label.set_fontsize(8)
        t.on_submit(accion)
        return ax, t

    def _boton(self, rect, texto, accion):
        """Botón que solo responde al clic izquierdo (el derecho abre el código)."""
        ax = self.fig.add_axes(rect)
        b = Button(ax, texto, color="#edf2f7", hovercolor="#cbd5e0")
        b.label.set_fontsize(8.3)

        def solo_izquierdo(evento):
            if getattr(evento, "button", 1) == 1:
                accion()
        b.on_clicked(solo_izquierdo)
        return ax, b

    def _crear_controles(self):
        c = Configuracion()
        self.ax_for, self.r_for = self._radio([0.01, 0.80, 0.11, 0.105], "Formación",
                                              list(pl.FORMACIONES), self.formacion, self._elegir_formacion)
        self.ax_sel, self.r_sel = self._radio([0.01, 0.62, 0.11, 0.105], "Selección",
                                              list(SELECCIONES), c.seleccion, self._cambio_en_vivo)
        self.ax_cru, self.r_cru = self._radio([0.01, 0.44, 0.11, 0.105], "Cruce",
                                              list(CRUCES), c.cruce, self._cambio_en_vivo)
        self.ax_mut, self.r_mut = self._radio([0.01, 0.26, 0.11, 0.105], "Mutación",
                                              list(MUTACIONES), c.mutacion, self._cambio_en_vivo)

        self.fig.text(0.135, 0.905, "Parámetros", fontsize=9, weight="bold")
        self.ax_pob, self.s_pob = self._slider(0.875, "Población", 4, 300, POBLACION_VENTANA, 1)
        self.ax_gen, self.s_gen = self._slider(0.845, "Generaciones máx.", 10, 1000, c.generaciones, 10)
        self.ax_pc, self.s_pc = self._slider(0.815, "Prob. cruce", 0, 1, c.prob_cruce, 0.05, "%.2f")
        self.ax_pm, self.s_pm = self._slider(0.785, "Prob. mutación", 0, 1, c.prob_mutacion, 0.01, "%.2f")
        self.ax_eli, self.s_eli = self._slider(0.755, "Elitismo (%)", 0, 50, c.elitismo_pct, 1)
        self.ax_k, self.s_k = self._slider(0.725, "Tamaño torneo k", 2, 10, c.k_torneo, 1)
        self.ax_lam, self.s_lam = self._slider(0.695, "Peso química λ", 0, 10, LAMBDA_VENTANA, 0.5, "%.1f")
        self.fig.text(0.135, 0.66, "Parada y escenario", fontsize=9, weight="bold")
        self.ax_sm, self.s_sm = self._slider(0.632, "Parar sin mejora", 0, 300, SIN_MEJORA_VENTANA, 10)
        self.ax_umb, self.t_umb = self._caja(0.595, "Puntaje objetivo",
                                             f"{OBJETIVO_CONOCIDO[self.formacion]:g}", self._cambio_en_vivo)
        self.ax_sem, self.t_sem = self._caja(0.560, "Semilla", "7", lambda _: self.reiniciar())
        self.ax_les, self.t_les = self._caja(0.525, "Lesionados (#)", "", self._cambio_escenario)

        for s in (self.s_pc, self.s_pm, self.s_eli, self.s_k, self.s_sm, self.s_gen):
            s.on_changed(self._cambio_en_vivo)
        self.s_pob.on_changed(lambda _: self.reiniciar())
        self.s_lam.on_changed(self._cambio_escenario)

        self.ax_run, self.b_run = self._boton([0.135, 0.468, 0.09, 0.036], "▶ Ejecutar", self._ejecutar_pausar)
        self.ax_paso, self.b_paso = self._boton([0.23, 0.468, 0.09, 0.036], "Un paso", self._un_paso)
        self.ax_rei, self.b_rei = self._boton([0.135, 0.426, 0.09, 0.036], "Reiniciar", self.reiniciar)
        self.ax_rep, self.b_rep = self._boton([0.23, 0.426, 0.09, 0.036], "Repetir evolución", self._repetir)
        # Tres ventanas aparte en la misma fila: la plantilla, los genes de cada mejora y el resumen.
        self.ax_pla, self.b_pla = self._boton([0.135, 0.384, 0.058, 0.036], "Plantilla", self._ver_plantilla)
        self.ax_vgen, self.b_vgen = self._boton([0.1985, 0.384, 0.058, 0.036], "Genes", self._ver_genes)
        self.ax_res, self.b_res = self._boton([0.262, 0.384, 0.058, 0.036], "Resumen", self._ver_resumen)

        self.ax_txt = self.fig.add_axes([0.135, 0.03, 0.2, 0.34])
        self.ax_txt.axis("off")
        self.texto = self.ax_txt.text(0, 1, "", va="top", ha="left", fontsize=7.4, family="monospace",
                                      transform=self.ax_txt.transAxes)
        # Titulares de la mejor alineación, debajo de «Mutación» (columna izquierda).
        self.ax_tit = self.fig.add_axes([0.01, 0.03, 0.115, 0.215])
        self.ax_tit.axis("off")
        self.titulares = self.ax_tit.text(0, 1, "", va="top", ha="left", fontsize=7.4, family="monospace",
                                          transform=self.ax_tit.transAxes)

    def _crear_graficas(self):
        self.ax_can = self.fig.add_axes([0.345, 0.095, 0.27, 0.835])
        self.ax_ver = self.fig.add_axes([0.40, 0.035, 0.18, 0.02])
        self.s_ver = Slider(self.ax_ver, "Ver generación", 0, 1, valinit=0, valstep=1, valfmt="%d")
        self.s_ver.label.set_fontsize(8)
        self.s_ver.valtext.set_fontsize(8)
        self.s_ver.on_changed(self._ver_generacion)
        self.ax_con = self.fig.add_axes([0.665, 0.715, 0.32, 0.215])
        self.ax_gen_pob = self.fig.add_axes([0.665, 0.415, 0.32, 0.21])
        self.ax_mat = self.fig.add_axes([0.71, 0.04, 0.275, 0.30])
        self._tabla_destinos()

    # ---------------------------------------------------------- problema y AG
    def _lesionados(self):
        numeros = []
        for parte in self.t_les.text.replace(";", ",").split(","):
            parte = parte.strip()
            if parte.isdigit() and 1 <= int(parte) <= pl.N_JUGADORES:
                numeros.append(int(parte) - 1)
        return numeros

    def _numero(self, caja, defecto):
        try:
            return float(caja.text.replace(",", "."))
        except ValueError:
            return defecto

    def _leer_configuracion(self):
        return Configuracion(
            tam_poblacion=int(self.s_pob.val), generaciones=int(self.s_gen.val),
            seleccion=self.r_sel.value_selected, k_torneo=int(self.s_k.val),
            cruce=self.r_cru.value_selected, prob_cruce=float(self.s_pc.val),
            mutacion=self.r_mut.value_selected, prob_mutacion=float(self.s_pm.val),
            elitismo_pct=float(self.s_eli.val), sin_mejora=int(self.s_sm.val),
            umbral=self._numero(self.t_umb, 0.0),
            semilla=abs(int(self._numero(self.t_sem, 7))))      # NumPy no acepta semillas negativas

    def reiniciar(self):
        self._pausar()
        self.problema = ProblemaAlineacion(self.formacion, float(self.s_lam.val), self._lesionados())
        self.referencia = optimo_sin_quimica(self.problema)
        self.ag = AlgoritmoGenetico(self.problema, self._leer_configuracion())
        self.gen_vista = None
        self._dibujar()

    # -------------------------------------------------------------- callbacks
    def _elegir_formacion(self, etiqueta):
        self.formacion = etiqueta
        self._ajustar_objetivo()
        self.reiniciar()

    def _cambio_escenario(self, _=None):
        """Cambió λ o los lesionados: el problema es otro y se empieza de nuevo."""
        self._ajustar_objetivo()
        self.reiniciar()

    def _ajustar_objetivo(self):
        """Pone como puntaje objetivo el mejor conocido de la formación.

        Solo se conoce para λ = 3 y sin lesionados; en otro caso se apaga (0),
        porque un objetivo equivocado haría parar antes de tiempo o nunca.
        """
        conocido = abs(float(self.s_lam.val) - LAMBDA_VENTANA) < 1e-9 and not self._lesionados()
        valor = f"{OBJETIVO_CONOCIDO[self.formacion]:g}" if conocido else "0"
        self.t_umb.eventson = False
        self.t_umb.set_val(valor)
        self.t_umb.eventson = True

    def _cambio_en_vivo(self, _=None):
        """Operadores, probabilidades, elitismo y parada cambian SIN reiniciar."""
        n = self._leer_configuracion()
        self.ag.actualizar_parametros(
            generaciones=n.generaciones, seleccion=n.seleccion, k_torneo=n.k_torneo, cruce=n.cruce,
            prob_cruce=n.prob_cruce, mutacion=n.mutacion, prob_mutacion=n.prob_mutacion,
            elitismo_pct=n.elitismo_pct, sin_mejora=n.sin_mejora, umbral=n.umbral)
        self._actualizar_texto()
        self.fig.canvas.draw_idle()

    def _ver_generacion(self, valor):
        if self.corriendo or self.repitiendo:
            return
        g = int(valor)
        self.gen_vista = None if g >= self.ag.generacion else g
        self._dibujar()

    # ---------------------------------------------------------------- correr
    def _ejecutar_pausar(self):
        if self.corriendo:
            self._pausar()
        else:
            self.repitiendo = False
            if self.ag.terminado():
                self.reiniciar()
            self.gen_vista = None
            self.corriendo = True
            self.b_run.label.set_text("❚❚ Pausa")
            self.timer.start()
        self.fig.canvas.draw_idle()

    def _pausar(self):
        self.corriendo = False
        self.repitiendo = False
        if hasattr(self, "timer"):
            self.timer.stop()
        if hasattr(self, "b_run"):
            self.b_run.label.set_text("▶ Ejecutar")

    def _repetir(self):
        """Reproduce la evolución guardada desde la generación 0."""
        self._pausar()
        self.repitiendo = True
        self.gen_vista = 0
        self.timer.start()

    def _tic(self):
        if self.repitiendo:
            self.gen_vista += max(1, self.ag.generacion // 60)
            if self.gen_vista >= self.ag.generacion:
                self._pausar()
                self.gen_vista = None
            self._dibujar()
            return
        if not self.corriendo:
            return
        if self.ag.terminado():
            self._pausar()
            self._dibujar()
            return
        self.ag.paso()
        self._dibujar()

    def _un_paso(self):
        self._pausar()
        if not self.ag.terminado():
            self.ag.paso()
        self.gen_vista = None
        self._dibujar()

    # ----------------------------------------------------------------- dibujo
    def _dibujar(self):
        ag, h = self.ag, self.ag.historial
        g = ag.generacion if self.gen_vista is None else self.gen_vista
        mejor = h.mejores[g]
        titulo = (f"Mejor alineación de la generación {g}  ·  puntaje {h.mejor[g]:.1f}")
        graficas.dibujar_cancha(self.ax_can, self.problema, mejor, titulo)
        graficas.dibujar_convergencia(self.ax_con, h, self.problema.evaluar(self.referencia["cromosoma"])[0],
                                      gen_marcada=g if self.gen_vista is not None else None)
        graficas.dibujar_genes(self.ax_gen_pob, self.problema, h.poblaciones[g], h.puntajes[g],
                               n_elite(ag.cfg), generacion=g)
        graficas.dibujar_matriz(self.ax_mat, self.problema, mejor)
        self.s_ver.eventson = False
        self.s_ver.valmax = max(1, ag.generacion)
        self.s_ver.ax.set_xlim(0, self.s_ver.valmax)
        self.s_ver.set_val(g)
        self.s_ver.eventson = True
        self._actualizar_texto()
        self.fig.canvas.draw_idle()

    def _texto(self):
        ag, cfg = self.ag, self.ag.cfg
        d = self.problema.detalle(ag.mejor)
        ref = self.problema.evaluar(self.referencia["cromosoma"])[0]
        # Compacto a propósito (12 líneas, hasta unos 50 caracteres): así cabe
        # con letra legible también en un portátil de 1366 × 768.
        parada = [f"{cfg.generaciones} gen."]
        if cfg.sin_mejora:
            parada.append(f"{cfg.sin_mejora} sin mejora")
        if cfg.umbral > 0:
            parada.append(f"objetivo {cfg.umbral:g}")
        lineas = [f"Generación {ag.generacion} / {cfg.generaciones}",
                  f"Paró: {ag.motivo_parada()}" if ag.terminado() else "",
                  f"{cfg.seleccion} · {cfg.cruce} · {cfg.mutacion}",
                  f"Población {cfg.tam_poblacion} · élite {n_elite(cfg)} · "
                  f"pc {cfg.prob_cruce:.2f} · pm {cfg.prob_mutacion:.2f}",
                  "Parada: " + " · ".join(parada),
                  f"Diversidad {ag.historial.diversidad[-1]:.1f} % de genes distintos",
                  "",
                  f"MEJOR PUNTAJE  {d['total']:8.1f}  (gen. {ag.gen_mejor})",
                  f"  rendimiento  {d['rendimiento']:8.1f}",
                  f"  química      {d['quimica']:+8.1f}",
                  f"  penalización {-d['penalizacion']:+8.1f}",
                  f"Húngaro + química {ref:7.1f}"]
        return "\n".join(lineas)

    def _texto_titulares(self):
        """Los 11 titulares de la mejor alineación: puesto, jugador, rendimiento y línea natural."""
        lineas = ["TITULARES", "(* = fuera de su línea)"]
        for f in self.problema.detalle(self.ag.mejor)["titulares"]:
            marca = "*" if graficas.fuera_de_posicion(pl.NOMBRES.index(f["jugador"]), f["rol"]) else " "
            lineas.append(f"{f['puesto']:4s}{marca}{f['corto']:9s}{f['rendimiento']:5.1f} {f['natural']}")
        return "\n".join(lineas)

    def _actualizar_texto(self):
        self.texto.set_text(self._texto())
        self.titulares.set_text(self._texto_titulares())
        self._ajustar_letras()

    def _ajustar_letras(self):
        self._ajustar_letra(self.texto, 0.345, 8.0)   # el panel llega hasta la cancha
        self._ajustar_letra(self.titulares, 0.135)    # los titulares, hasta el panel

    def _ajustar_letra(self, texto, borde_derecho, maxima=7.4):
        """Pone en `texto` la letra más grande (hasta `maxima`) que quepa.

        El texto va desde su borde superior hasta el pie de la ventana y, a lo
        ancho, hasta `borde_derecho` (fracción de la figura). Con la lupa de
        ajustar_a_pantalla las letras crecen y el texto podría salirse por
        abajo o invadir lo de al lado; aquí se reduce lo justo.
        """
        # El ancho de la letra monoespaciada se redondea a píxeles enteros, así
        # que el tamaño del texto crece a saltos: se mide de verdad y se baja
        # un 4 % cada vez hasta que quepa (casi siempre cabe a la primera).
        letra = maxima
        while letra > 4:
            texto.set_fontsize(letra)
            try:
                caja = texto.get_window_extent(self.fig.canvas.get_renderer())
            except Exception:
                return                               # backend sin renderer: se deja como está
            ancho_px, alto_px = self.fig.bbox.width, self.fig.bbox.height
            if caja.y0 >= 0.004 * alto_px and caja.x1 <= (borde_derecho - 0.004) * ancho_px:
                return
            letra *= 0.96

    def _ver_plantilla(self):
        fig = graficas.fig_plantilla()
        fig.canvas.mpl_connect("button_press_event", lambda e: self._abrir_si_derecho(
            e, "plantilla.py", "JUGADORES"))
        ajustar_a_pantalla(fig, 11, 9.4)
        fig.show()

    def _ventana_con_guardar(self, fig, nombre_base):
        """Pone en una ventana aparte el botón «Guardar imagen» (guarda la principal y esa ventana)."""
        ax_b = fig.add_axes([0.87, 0.012, 0.115, 0.045])
        boton = Button(ax_b, "Guardar imagen", color="#edf2f7", hovercolor="#e2e8f0")
        boton.label.set_fontsize(8.5)
        aviso = fig.text(0.40, 0.025, "", fontsize=8, color=graficas.AZUL)

        def guardar(_):
            principal = self._guardar()
            nombre = f"{nombre_base}_generacion_{self.ag.generacion}.png"
            fig.savefig(abrir_codigo.ruta(nombre), dpi=110)
            aviso.set_text(f"Guardadas: {principal} y {nombre}")
            fig.canvas.draw_idle()
        boton.on_clicked(guardar)
        return boton

    def _ver_genes(self):
        """Ventana aparte con el mejor cromosoma de cada generación en que mejoró la alineación.

        El deslizador elige hasta qué generación se muestran las mejoras.
        """
        h, ag = self.ag.historial, self.ag
        fig = plt.figure(figsize=(10, 8.6))
        ax_mej = fig.add_axes([0.16, 0.14, 0.81, 0.79])
        ax_sl = fig.add_axes([0.16, 0.035, 0.45, 0.03])
        actual = self.gen_vista if self.gen_vista is not None else ag.generacion
        deslizador = Slider(ax_sl, "Hasta la generación", 0, max(1, ag.generacion), valinit=actual,
                            valstep=1, valfmt="%d")

        def dibujar(g):
            graficas.dibujar_mejores(ax_mej, self.problema, h, int(min(g, ag.generacion)))
            fig.canvas.draw_idle()

        deslizador.on_changed(dibujar)
        boton = self._ventana_con_guardar(fig, "genes")
        fig.canvas.mpl_connect("button_press_event", lambda e: self._abrir_si_derecho(
            e, "graficas.py", "dibujar_mejores"))
        self._ventana_genes = (fig, deslizador, boton)   # sin referencia los controles dejan de responder
        ajustar_a_pantalla(fig, 10, 8.6)
        dibujar(actual)
        fig.show()
        return fig

    def _ver_resumen(self):
        """Ventana aparte con la evolución, la diversidad, las mejoras y los datos de la corrida."""
        fig = graficas.fig_resumen(self.ag, self.problema.evaluar(self.referencia["cromosoma"])[0])
        boton = self._ventana_con_guardar(fig, "resumen")
        fig.canvas.mpl_connect("button_press_event", lambda e: self._abrir_si_derecho(
            e, "graficas.py", "fig_resumen"))
        self._resumen = (fig, boton)
        ajustar_a_pantalla(fig, 14, 8.6)
        fig.show()
        return fig

    def _guardar(self, nombre=None):
        nombre = nombre or f"captura_generacion_{self.ag.generacion}.png"
        self.fig.savefig(abrir_codigo.ruta(nombre), dpi=110)
        self.estado.set_text(f"Imagen guardada: {nombre}")
        self.fig.canvas.draw_idle()
        return nombre

    # ----------------------------------------------------------- clic derecho
    def _tabla_destinos(self):
        """Control -> (archivo, función principal, [otras funciones que intervienen])."""
        self.destinos = {
            self.ax_for: ("plantilla.py", "FORMACIONES", [("plantilla.py", "vecinos")]),
            self.ax_pob: ("algoritmo.py", "AlgoritmoGenetico.reiniciar", []),
            self.ax_gen: ("algoritmo.py", "AlgoritmoGenetico.terminado", []),
            self.ax_pc: ("operadores.py", "cruzar_poblacion", []),
            self.ax_pm: ("operadores.py", "mutar_poblacion", []),
            self.ax_eli: ("algoritmo.py", "n_elite", [("algoritmo.py", "AlgoritmoGenetico.paso")]),
            self.ax_k: ("operadores.py", "seleccion_torneo", []),
            self.ax_lam: ("evaluacion.py", "ProblemaAlineacion.puntaje_quimica", [("plantilla.py", "matriz_quimica")]),
            self.ax_sm: ("algoritmo.py", "AlgoritmoGenetico.terminado", []),
            self.ax_umb: ("algoritmo.py", "AlgoritmoGenetico.terminado", []),
            self.ax_sem: ("algoritmo.py", "AlgoritmoGenetico.reiniciar", []),
            self.ax_les: ("evaluacion.py", "ProblemaAlineacion.penalizacion", []),
            self.ax_run: ("algoritmo.py", "AlgoritmoGenetico.paso", [("interfaz.py", "Simulador._tic")]),
            self.ax_paso: ("algoritmo.py", "AlgoritmoGenetico.paso", []),
            self.ax_rei: ("algoritmo.py", "AlgoritmoGenetico.reiniciar", []),
            self.ax_rep: ("interfaz.py", "Simulador._repetir", []),
            self.ax_pla: ("plantilla.py", "JUGADORES", [("plantilla.py", "ROLES")]),
            self.ax_vgen: ("graficas.py", "dibujar_mejores", [("interfaz.py", "Simulador._ver_genes")]),
            self.ax_res: ("graficas.py", "fig_resumen", [("graficas.py", "cambios_entre")]),
            self.ax_txt: ("evaluacion.py", "ProblemaAlineacion.detalle", [("interfaz.py", "Simulador._texto")]),
            self.ax_tit: ("evaluacion.py", "ProblemaAlineacion.detalle",
                          [("interfaz.py", "Simulador._texto_titulares")]),
            self.ax_can: ("graficas.py", "dibujar_cancha", [("evaluacion.py", "ProblemaAlineacion.evaluar")]),
            self.ax_ver: ("interfaz.py", "Simulador._ver_generacion", []),
            self.ax_con: ("graficas.py", "dibujar_convergencia", [("evaluacion.py", "optimo_sin_quimica")]),
            self.ax_gen_pob: ("graficas.py", "dibujar_genes", [("algoritmo.py", "AlgoritmoGenetico.diversidad")]),
            self.ax_mat: ("evaluacion.py", "matriz_rendimiento", [("evaluacion.py", "rendimiento_base")]),
        }

    def _opcion_bajo_clic(self, radio, evento):
        for etiqueta in radio.labels:
            caja = etiqueta.get_window_extent()
            if caja.y0 - 3 <= evento.y <= caja.y1 + 3:
                return etiqueta.get_text()
        return None

    def destino_de(self, evento):
        ax = evento.inaxes
        por_opcion = {
            self.ax_sel: (self.r_sel, SELECCIONES),
            self.ax_cru: (self.r_cru, CRUCES),
            self.ax_mut: (self.r_mut, MUTACIONES),
        }
        if ax in por_opcion:
            radio, tabla = por_opcion[ax]
            opcion = self._opcion_bajo_clic(radio, evento)
            if opcion is not None:
                return "operadores.py", tabla[opcion].__name__, [("algoritmo.py", "AlgoritmoGenetico.paso")]
            return "interfaz.py", "Simulador._leer_configuracion", []
        if ax in self.destinos:
            return self.destinos[ax]
        return None

    def _abrir_si_derecho(self, evento, archivo, nombre, otros=()):
        if evento.button != 3:
            return
        abierto, donde = abrir_codigo.abrir(archivo, nombre)
        self.estado.set_text(f"Abriendo {donde} · {nombre}" if abierto
                             else f"No se pudo abrir VS Code. Ve a {donde} ({nombre})")
        self.fig.canvas.draw_idle()

    def _clic_derecho(self, evento):
        if evento.button != 3 or evento.inaxes is None:
            return
        barra = getattr(self.fig.canvas, "toolbar", None)
        if barra is not None and getattr(barra, "mode", ""):
            return
        destino = self.destino_de(evento)
        if destino is None:
            return
        archivo, nombre, otros = destino
        if getattr(evento, "key", None) == "shift" and otros:
            (archivo, nombre), otros = otros[0], []
        abierto, donde = abrir_codigo.abrir(archivo, nombre)
        mensaje = f"Abriendo {donde} · {nombre}" if abierto else f"No se pudo abrir VS Code. Ve a {donde} ({nombre})"
        if otros:
            mensaje += "   También: " + ", ".join(f"{a} › {n}" for a, n in otros)
        self.estado.set_text(mensaje)
        self.ultimo_destino = (archivo, nombre)
        self.fig.canvas.draw_idle()


# Cuánto se agrandan las letras (y las líneas) en una pantalla pequeña respecto
# a encoger la ventana de forma proporcional. Con más de 1.25 los títulos y los
# controles empiezan a chocar en 1366 × 768.
AUMENTO_LETRA = 1.25


def ajustar_a_pantalla(fig, ancho_pulg=16, alto_pulg=9, aumento=AUMENTO_LETRA):
    """Acomoda la ventana si la pantalla es más pequeña que la figura.

    Encoger todo en la misma proporción deja letras de 6 o 7 píxeles en un
    portátil de 1366 × 768. Por eso, en una pantalla pequeña:
      * se quita la barra de herramientas de Matplotlib (la ventana ya tiene
        «Guardar imagen» dentro de «Genes» y «Resumen») y se maximiza la ventana, para usar toda la pantalla;
      * el dpi baja menos de lo que haría falta para encoger (lupa = aumento):
        las letras y las líneas, que se miden en puntos, quedan más grandes,
        y las posiciones de los controles, que son fracciones de la figura,
        no cambian.
    Devuelve la escala de la pantalla (1 = sin cambios).
    """
    try:
        ventana = fig.canvas.manager.window
        ancho = ventana.winfo_screenwidth()
        alto = ventana.winfo_screenheight()
    except Exception:
        return 1.0                                   # otro backend o sin pantalla: no se toca
    dpi = fig.dpi
    escala = min(1.0, (ancho - 40) / (ancho_pulg * dpi), (alto - 120) / (alto_pulg * dpi))
    if escala < 1.0:
        factor = escala * min(aumento, 1 / escala)   # nunca más grande que en pantalla grande
        # _original_dpi es el dpi sin el factor de pantallas de alta densidad;
        # Matplotlib lo vuelve a multiplicar si cambia ese factor.
        if hasattr(fig, "_original_dpi"):
            fig._original_dpi *= factor
        fig.set_dpi(dpi * factor)
        barra = getattr(fig.canvas, "toolbar", None)
        if barra is not None:
            try:
                barra.pack_forget()
            except Exception:
                pass
        # Tamaño de la figura = espacio libre de la pantalla (sin barra de tareas ni título).
        fig.set_size_inches((ancho - 16) / fig.dpi, (alto - 80) / fig.dpi, forward=True)
        try:
            ventana.state("zoomed")                  # Windows: maximizada
        except Exception:
            try:
                ventana.geometry("+0+0")             # otros sistemas: esquina superior izquierda
            except Exception:
                pass
    return escala


def lanzar():
    sim = Simulador()
    ajustar_a_pantalla(sim.fig)
    plt.show()
    return sim
