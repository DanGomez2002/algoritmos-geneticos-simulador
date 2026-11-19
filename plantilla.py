"""
plantilla.py
------------
Los datos del problema: los 18 jugadores de la plantilla con sus
estadísticas, el perfil de cada rol (qué atributos pesan en cada puesto),
las formaciones y la química entre jugadores.

Todo es inventado con criterio de entrenador; los valores se pueden cambiar
aquí sin tocar el resto del programa.

Escala de los atributos: 1 a 99, como en los videojuegos de fútbol.
"""

import numpy as np

ATRIBUTOS = ["Ritmo", "Tiro", "Pase", "Regate", "Defensa", "Físico", "Portería"]

# ---------------------------------------------------------------------------
# JUGADORES
# nombre, corto (para la cancha), posición natural, pie (D, I o A = ambos),
# ritmo, tiro, pase, regate, defensa, físico, portería, estado físico (%)
# ---------------------------------------------------------------------------
JUGADORES = [
    ("Mateo Ríos",          "Ríos",      "POR", "D", 45, 20, 55, 30, 25, 70, 82, 100),
    ("Samuel Cuéllar",      "Cuéllar",   "POR", "I", 50, 18, 60, 35, 22, 65, 76, 100),
    ("Jhon Fredy Perdomo",  "Perdomo",   "DFC", "D", 60, 35, 58, 45, 84, 82, 12, 100),
    ("Kevin Andrés Silva",  "Silva",     "DFC", "I", 68, 30, 62, 50, 80, 76, 10, 95),
    ("Brayan Castaño",      "Castaño",   "DFC", "D", 62, 45, 70, 55, 76, 80, 11, 100),
    ("Cristian Vargas",     "Vargas",    "LD",  "D", 82, 40, 68, 70, 74, 66, 10, 100),
    ("Luis Ángel Murcia",   "Murcia",    "LI",  "I", 80, 42, 70, 72, 72, 62, 9, 90),
    ("Óscar Trujillo",      "Trujillo",  "LD",  "A", 74, 38, 64, 62, 70, 70, 10, 100),
    ("Daniel Ospina",       "Ospina",    "MCD", "D", 60, 55, 76, 62, 78, 82, 12, 100),
    ("Santiago Restrepo",   "Restrepo",  "MC",  "D", 66, 66, 86, 80, 60, 68, 10, 100),
    ("Felipe Cardona",      "Cardona",   "MC",  "I", 70, 70, 80, 78, 62, 70, 11, 95),
    ("Nicolás Zapata",      "Zapata",    "MCO", "I", 74, 80, 82, 86, 40, 58, 9, 75),
    ("Esteban Monroy",      "Monroy",    "ED",  "D", 88, 70, 72, 82, 45, 60, 8, 100),
    ("Alexis Guzmán",       "Guzmán",    "EI",  "I", 86, 72, 70, 84, 42, 58, 9, 100),
    ("Yeison Mosquera",     "Mosquera",  "DC",  "D", 78, 88, 60, 74, 35, 84, 10, 100),
    ("Julián Patiño",       "Patiño",    "DC",  "D", 72, 76, 66, 68, 72, 85, 11, 100),
    ("Diego Bermúdez",      "Bermúdez",  "ED",  "A", 90, 74, 64, 80, 30, 55, 8, 85),
    ("Harold Calderón",     "Calderón",  "MC",  "I", 72, 62, 74, 70, 70, 74, 12, 100),
]

NOMBRES = [j[0] for j in JUGADORES]
CORTOS = [j[1] for j in JUGADORES]
NATURAL = [j[2] for j in JUGADORES]
PIE = [j[3] for j in JUGADORES]
STATS = np.array([j[4:11] for j in JUGADORES], float)          # (18 x 7)
ESTADO = np.array([j[11] for j in JUGADORES], float)            # (18,)
N_JUGADORES = len(JUGADORES)
ES_PORTERO = np.array([n == "POR" for n in NATURAL])

# ---------------------------------------------------------------------------
# ROLES: peso de cada atributo en cada rol (cada fila suma 1).
# Ritmo, Tiro, Pase, Regate, Defensa, Físico, Portería
# ---------------------------------------------------------------------------
ROLES = {
    "POR": [0.05, 0.00, 0.05, 0.00, 0.00, 0.05, 0.85],   # casi todo es atajar
    "DFC": [0.10, 0.00, 0.10, 0.05, 0.50, 0.25, 0.00],   # marca y juego aéreo
    "LAT": [0.25, 0.00, 0.20, 0.10, 0.35, 0.10, 0.00],   # defender y subir por la banda
    "CAR": [0.30, 0.00, 0.20, 0.15, 0.25, 0.10, 0.00],   # carrilero: más ida y vuelta
    "MCD": [0.05, 0.05, 0.25, 0.05, 0.35, 0.25, 0.00],   # recuperar y entregar
    "MC":  [0.10, 0.10, 0.35, 0.15, 0.15, 0.15, 0.00],   # distribuir el juego
    "MED": [0.30, 0.05, 0.25, 0.20, 0.10, 0.10, 0.00],   # volante por banda
    "EXT": [0.30, 0.20, 0.15, 0.30, 0.00, 0.05, 0.00],   # desborde y gol
    "DC":  [0.15, 0.45, 0.05, 0.15, 0.00, 0.20, 0.00],   # definición y choque
}
ROL_NOMBRE = {
    "POR": "Portero", "DFC": "Defensa central", "LAT": "Lateral", "CAR": "Carrilero",
    "MCD": "Volante de marca", "MC": "Volante central", "MED": "Volante por banda",
    "EXT": "Extremo", "DC": "Delantero centro",
}

# ---------------------------------------------------------------------------
# FORMACIONES: lista de puestos (sigla, rol, lado, x, y).
# x: 0 = banda izquierda, 100 = banda derecha; y: 0 = arco propio, 100 = rival.
# lado: "D", "I" o "" (centro). El orden de esta lista es el orden de los
# primeros 11 genes del cromosoma.
# ---------------------------------------------------------------------------
FORMACIONES = {
    "4-4-2": [
        ("POR", "POR", "", 50, 4),
        ("LD", "LAT", "D", 86, 28), ("DFC", "DFC", "", 62, 23),
        ("DFC", "DFC", "", 38, 23), ("LI", "LAT", "I", 14, 28),
        ("MD", "MED", "D", 86, 54), ("MC", "MC", "", 62, 48),
        ("MC", "MC", "", 38, 48), ("MI", "MED", "I", 14, 54),
        ("DC", "DC", "", 62, 80), ("DC", "DC", "", 38, 80),
    ],
    "4-3-3": [
        ("POR", "POR", "", 50, 4),
        ("LD", "LAT", "D", 86, 28), ("DFC", "DFC", "", 62, 23),
        ("DFC", "DFC", "", 38, 23), ("LI", "LAT", "I", 14, 28),
        ("MCD", "MCD", "", 50, 40), ("MC", "MC", "", 70, 54),
        ("MC", "MC", "", 30, 54), ("ED", "EXT", "D", 84, 78),
        ("DC", "DC", "", 50, 84), ("EI", "EXT", "I", 16, 78),
    ],
    "3-5-2": [
        ("POR", "POR", "", 50, 4),
        ("DFC", "DFC", "", 72, 25), ("DFC", "DFC", "", 50, 22),
        ("DFC", "DFC", "", 28, 25), ("CD", "CAR", "D", 88, 50),
        ("MCD", "MCD", "", 50, 38), ("MC", "MC", "", 68, 56),
        ("MC", "MC", "", 32, 56), ("CI", "CAR", "I", 12, 50),
        ("DC", "DC", "", 62, 82), ("DC", "DC", "", 38, 82),
    ],
}
N_TITULARES = 11
DISTANCIA_VECINOS = 35          # dos puestos más cerca que esto son vecinos


def vecinos(formacion):
    """Pares de puestos (i, j) que juegan cerca: ahí cuenta la química."""
    p = FORMACIONES[formacion]
    pares = []
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            d = np.hypot(p[i][3] - p[j][3], p[i][4] - p[j][4])
            if d < DISTANCIA_VECINOS:
                pares.append((i, j))
    return pares


# ---------------------------------------------------------------------------
# QUÍMICA entre pares de jugadores: +2 se entienden muy bien, +1 bien,
# -1 regular, -2 no se entienden. Solo suma si juegan en puestos vecinos.
# ---------------------------------------------------------------------------
_PARES_QUIMICA = [
    ("Perdomo", "Silva", +2),      # dupla de centrales desde las divisiones menores
    ("Vargas", "Monroy", +2),      # lateral y extremo derecho que se buscan de memoria
    ("Murcia", "Guzmán", +2),      # lo mismo por la izquierda
    ("Restrepo", "Zapata", +2),    # el que distribuye y el que recibe entre líneas
    ("Ospina", "Restrepo", +1),
    ("Mosquera", "Patiño", +1),
    ("Ríos", "Perdomo", +1),       # el portero ordena bien a su central de confianza
    ("Cardona", "Zapata", -2),     # se pelean la pelota y los tiros libres
    ("Castaño", "Silva", -1),
    ("Bermúdez", "Vargas", -1),
]


def matriz_quimica():
    Q = np.zeros((N_JUGADORES, N_JUGADORES))
    for a, b, v in _PARES_QUIMICA:
        i, j = CORTOS.index(a), CORTOS.index(b)
        Q[i, j] = Q[j, i] = v
    return Q


QUIMICA = matriz_quimica()
