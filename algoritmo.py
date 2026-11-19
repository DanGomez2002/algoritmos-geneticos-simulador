"""
algoritmo.py
------------
El ciclo del algoritmo genético para la alineación.

    1. Primera generación: N alineaciones aleatorias (permutaciones al azar).
    2. Evaluar el puntaje de cada cromosoma (evaluacion.py).
    3. Elitismo: el porcentaje de los mejores pasa INTACTO a la siguiente
       generación (garantiza que el mejor puntaje nunca baje).
    4. Selección de padres: los de mejor puntaje tienen más probabilidad.
    5. Cruce por parejas con probabilidad pc (PMX, OX o ciclos).
    6. Mutación de cada hijo con probabilidad pm (porcentaje pequeño).
    7. La nueva generación reemplaza a la anterior. Se repite hasta una
       condición de parada: número máximo de generaciones, generaciones sin
       mejora o alcanzar un puntaje objetivo (umbral).
"""

from dataclasses import dataclass, field, replace

import numpy as np

import plantilla as pl
from operadores import SELECCIONES, cruzar_poblacion, mutar_poblacion


@dataclass
class Configuracion:
    tam_poblacion: int = 50
    generaciones: int = 200
    seleccion: str = "Ruleta"
    k_torneo: int = 3
    cruce: str = "PMX"
    prob_cruce: float = 0.9
    mutacion: str = "Intercambio"
    prob_mutacion: float = 0.10
    elitismo_pct: float = 10.0          # % de la población que pasa intacta
    sin_mejora: int = 0                 # 0 = no parar por estancamiento
    umbral: float = 0.0                 # 0 = sin puntaje objetivo
    semilla: int = 7


def n_elite(cfg):
    """Cuántos cromosomas pasan por elitismo (al menos 1 si el porcentaje es > 0)."""
    if cfg.elitismo_pct <= 0:
        return 0
    return int(min(cfg.tam_poblacion - 1, max(1, round(cfg.tam_poblacion * cfg.elitismo_pct / 100))))


def aptitud_desde_puntaje(puntajes):
    """Lleva los puntajes a valores >= 0 para la ruleta (escalamiento por ventana).

        a_i = (p_i − p_peor) + 0,01 · (p_mejor − p_peor)

    El peor queda con una probabilidad pequeña pero no nula de ser padre.
    """
    p = np.asarray(puntajes, float)
    rango = p.max() - p.min()
    if rango <= 1e-12:
        return np.ones_like(p)
    return (p - p.min()) + 0.01 * rango


@dataclass
class Historial:
    mejor: list = field(default_factory=list)          # mejor puntaje de cada generación
    media: list = field(default_factory=list)
    peor: list = field(default_factory=list)
    mejor_global: list = field(default_factory=list)
    diversidad: list = field(default_factory=list)     # genes distintos del mejor, en promedio (%)
    poblaciones: list = field(default_factory=list)    # cromosomas de cada generación, ordenados
    puntajes: list = field(default_factory=list)
    mejores: list = field(default_factory=list)        # mejor cromosoma de cada generación


class AlgoritmoGenetico:

    def __init__(self, problema, config=None, guardar_poblaciones=True):
        self.problema = problema
        self.cfg = config or Configuracion()
        self.guardar = guardar_poblaciones
        self.reiniciar()

    def reiniciar(self):
        self.rng = np.random.default_rng(self.cfg.semilla)
        N, L = self.cfg.tam_poblacion, pl.N_JUGADORES
        self.poblacion = np.array([self.rng.permutation(L) for _ in range(N)])
        self.generacion = 0
        self.evaluaciones = 0
        self.historial = Historial()
        self.mejor_puntaje = -np.inf
        self.mejor = None
        self.gen_mejor = 0
        self._sin_mejora = 0
        self._evaluar()

    def actualizar_parametros(self, **cambios):
        """Cambia parámetros sin reiniciar (cruce, mutación, selección, pc, pm, elitismo, parada)."""
        self.cfg = replace(self.cfg, **cambios)

    def _evaluar(self):
        self.puntajes = self.problema.evaluar(self.poblacion)
        self.evaluaciones += self.puntajes.size
        orden = np.argsort(-self.puntajes)                    # del mejor al peor
        self.poblacion = self.poblacion[orden]
        self.puntajes = self.puntajes[orden]
        if self.puntajes[0] > self.mejor_puntaje + 1e-9:
            self.mejor_puntaje = float(self.puntajes[0])
            self.mejor = self.poblacion[0].copy()
            self.gen_mejor = self.generacion
            self._sin_mejora = 0
        else:
            self._sin_mejora += 1
        h = self.historial
        h.mejor.append(float(self.puntajes[0]))
        h.media.append(float(self.puntajes.mean()))
        h.peor.append(float(self.puntajes[-1]))
        h.mejor_global.append(self.mejor_puntaje)
        h.diversidad.append(self.diversidad())
        h.mejores.append(self.poblacion[0].copy())
        if self.guardar:
            h.poblaciones.append(self.poblacion.copy())
            h.puntajes.append(self.puntajes.copy())

    def diversidad(self):
        """Porcentaje promedio de genes (puestos) en que cada cromosoma difiere del mejor.

        100 % = nadie se parece al mejor; 0 % = toda la población es la misma alineación.
        """
        return float(100 * np.mean(self.poblacion != self.poblacion[0]))

    def paso(self):
        """Una generación completa."""
        cfg = self.cfg
        N = cfg.tam_poblacion
        e = n_elite(cfg)
        elite = self.poblacion[:e].copy()                     # ya está ordenada
        aptitud = aptitud_desde_puntaje(self.puntajes)
        idx = SELECCIONES[cfg.seleccion](aptitud, N - e, self.rng, k=cfg.k_torneo)
        padres = self.poblacion[idx]
        hijos = cruzar_poblacion(padres, cfg.prob_cruce, cfg.cruce, self.rng)
        hijos = mutar_poblacion(hijos, cfg.prob_mutacion, cfg.mutacion, self.rng)
        self.poblacion = np.vstack([elite, hijos]) if e else hijos
        self.generacion += 1
        self._evaluar()

    def terminado(self):
        cfg = self.cfg
        if self.generacion >= cfg.generaciones:
            return True
        if cfg.sin_mejora > 0 and self._sin_mejora >= cfg.sin_mejora:
            return True
        return cfg.umbral > 0 and self.mejor_puntaje >= cfg.umbral

    def motivo_parada(self):
        cfg = self.cfg
        if cfg.umbral > 0 and self.mejor_puntaje >= cfg.umbral:
            return "alcanzó el puntaje objetivo"
        if cfg.sin_mejora > 0 and self._sin_mejora >= cfg.sin_mejora:
            return f"{cfg.sin_mejora} generaciones sin mejorar"
        if self.generacion >= cfg.generaciones:
            return "llegó al máximo de generaciones"
        return ""

    def ejecutar(self, al_paso=None):
        while not self.terminado():
            self.paso()
            if al_paso:
                al_paso(self)
        return self.resultado()

    def resultado(self):
        return {"mejor_puntaje": self.mejor_puntaje, "mejor": self.mejor.copy(),
                "generaciones": self.generacion, "gen_mejor": self.gen_mejor,
                "evaluaciones": self.evaluaciones, "diversidad_final": self.historial.diversidad[-1]}
