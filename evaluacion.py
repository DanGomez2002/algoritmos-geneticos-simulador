"""
evaluacion.py
-------------
La función de evaluación (fitness): qué tan buena es una alineación.

Cromosoma = permutación de los 18 jugadores. Los genes 0 a 10 son los
puestos de la formación (en el orden de plantilla.FORMACIONES) y los genes
11 a 17 son el banco.

    puntaje = Σ rendimiento(jugador, puesto)            (los 11 titulares)
            + λ · Σ química(jugador_a, jugador_b)        (pares de puestos vecinos)
            − penalización por banco sin portero suplente (sano)
            − penalización por cada lesionado en la cancha

El rendimiento de un jugador en un puesto sale de sus estadísticas:

    rendimiento = (Σ_k peso_rol,k · estadística_k) · factor_pie · factor_físico

* factor_pie: en un puesto de banda, jugar con la pierna contraria
  (un diestro de lateral izquierdo) rinde 8 % menos; el ambidiestro no pierde.
* factor_físico = 0,7 + 0,3 · estado/100: un jugador al 75 % rinde 92,5 %.
"""

import numpy as np

import plantilla as pl

PENAL_PORTERO_BANCO = 40.0      # sin portero suplente en el banco
PENAL_LESIONADO = 300.0         # por cada lesionado que se ponga de titular
CASTIGO_PIE = 0.92


def rendimiento_base(jugador, rol):
    """Suma ponderada de las estadísticas del jugador según el rol."""
    return float(pl.STATS[jugador] @ np.array(pl.ROLES[rol]))


def factor_pie(jugador, lado):
    pie = pl.PIE[jugador]
    if lado == "" or pie == "A" or pie == lado:
        return 1.0
    return CASTIGO_PIE


def factor_fisico(jugador):
    return 0.7 + 0.3 * pl.ESTADO[jugador] / 100


def matriz_rendimiento(formacion):
    """Matriz R (18 jugadores × 11 puestos) con el rendimiento de cada uno en cada puesto."""
    puestos = pl.FORMACIONES[formacion]
    R = np.zeros((pl.N_JUGADORES, len(puestos)))
    for j in range(pl.N_JUGADORES):
        for p, (_, rol, lado, _, _) in enumerate(puestos):
            R[j, p] = rendimiento_base(j, rol) * factor_pie(j, lado) * factor_fisico(j)
    return R


class ProblemaAlineacion:
    """Todo lo necesario para evaluar alineaciones de una formación dada."""

    def __init__(self, formacion="4-4-2", lambda_quimica=3.0, lesionados=()):
        self.formacion = formacion
        self.puestos = pl.FORMACIONES[formacion]
        self.R = matriz_rendimiento(formacion)
        self.pares = np.array(pl.vecinos(formacion))
        self.lambda_quimica = float(lambda_quimica)
        self.lesionados = np.array(sorted(set(lesionados)), dtype=int)
        self.n = pl.N_JUGADORES

    # -- piezas del puntaje (vectorizadas: P es una matriz N x 18) ------------
    def puntaje_rendimiento(self, P):
        T = P[:, :pl.N_TITULARES]
        return self.R[T, np.arange(pl.N_TITULARES)].sum(axis=1)

    def puntaje_quimica(self, P):
        a = P[:, self.pares[:, 0]]
        b = P[:, self.pares[:, 1]]
        return self.lambda_quimica * pl.QUIMICA[a, b].sum(axis=1)

    def penalizacion(self, P):
        banco = P[:, pl.N_TITULARES:]
        # El suplente tiene que ser un portero sano: un lesionado en el banco no cuenta.
        portero_sano = pl.ES_PORTERO[banco] & ~np.isin(banco, self.lesionados)
        sin_portero = ~portero_sano.any(axis=1)
        pen = PENAL_PORTERO_BANCO * sin_portero
        if self.lesionados.size:
            en_cancha = np.isin(P[:, :pl.N_TITULARES], self.lesionados).sum(axis=1)
            pen = pen + PENAL_LESIONADO * en_cancha
        return pen

    def evaluar(self, P):
        """Puntaje de cada cromosoma (más alto = mejor alineación)."""
        P = np.atleast_2d(P)
        return self.puntaje_rendimiento(P) + self.puntaje_quimica(P) - self.penalizacion(P)

    def detalle(self, cromosoma):
        """Desglose de un cromosoma para mostrar en pantalla y en el informe."""
        c = np.asarray(cromosoma)
        P = c[None, :]
        filas = []
        for p, (sigla, rol, lado, _, _) in enumerate(self.puestos):
            j = int(c[p])
            filas.append({"puesto": sigla, "rol": rol, "jugador": pl.NOMBRES[j],
                          "corto": pl.CORTOS[j], "natural": pl.NATURAL[j],
                          "rendimiento": float(self.R[j, p])})
        quim = [(int(i), int(k), float(pl.QUIMICA[c[i], c[k]])) for i, k in self.pares
                if pl.QUIMICA[c[i], c[k]] != 0]
        return {
            "titulares": filas,
            "banco": [pl.NOMBRES[j] for j in c[pl.N_TITULARES:]],
            "rendimiento": float(self.puntaje_rendimiento(P)[0]),
            "quimica": float(self.puntaje_quimica(P)[0]),
            "penalizacion": float(self.penalizacion(P)[0]),
            "total": float(self.evaluar(P)[0]),
            "pares_quimica": quim,
        }


def optimo_sin_quimica(problema):
    """Referencia exacta con el método húngaro (scipy.optimize.linear_sum_assignment).

    Sin química la puntuación es una suma de valores jugador-puesto
    independientes: es el clásico problema de asignación y se resuelve de
    forma exacta. Con química deja de serlo (depende de pares de jugadores)
    y ahí es donde el algoritmo genético hace falta.
    """
    from scipy.optimize import linear_sum_assignment
    R = problema.R.copy()
    if problema.lesionados.size:
        R[problema.lesionados] = -1e6
    filas, cols = linear_sum_assignment(-R)
    cromosoma = np.empty(pl.N_TITULARES, int)
    cromosoma[cols] = filas
    resto = [j for j in range(pl.N_JUGADORES) if j not in set(filas)]
    completo = np.concatenate([cromosoma, resto])
    return {"cromosoma": completo, "rendimiento": float(R[filas, cols].sum())}
