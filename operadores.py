"""
operadores.py
-------------
Selección, cruce y mutación para cromosomas de PERMUTACIÓN.

En una permutación cada jugador aparece exactamente una vez. Por eso no
sirve el cruce de un punto de la codificación binaria: el hijo saldría con
jugadores repetidos y otros faltantes. Los cruces de aquí (PMX, OX y
ciclos) toman información de los dos padres SIN repetir a nadie, y las
mutaciones solo cambian el orden (nunca crean ni borran jugadores).

`aptitud`: vector donde más alto = mejor (siempre >= 0).
`rng`: generador de NumPy, para poder repetir una corrida con la misma semilla.
"""

import numpy as np


# ---------------------------------------------------------------------------
# SELECCIÓN DE PADRES (devuelven los índices de los elegidos)
# ---------------------------------------------------------------------------

def seleccion_ruleta(aptitud, n, rng, **_):
    """Ruleta: probabilidad de ser padre proporcional al puntaje, p_i = a_i / Σ a.

    Los mejores se cruzan más, pero todos tienen alguna probabilidad.
    """
    aptitud = np.asarray(aptitud, float)
    total = aptitud.sum()
    if total <= 0 or not np.isfinite(total):
        return rng.integers(0, aptitud.size, n)
    acumulada = np.cumsum(aptitud / total)
    acumulada[-1] = 1.0
    return np.searchsorted(acumulada, rng.random(n), side="right")


def seleccion_torneo(aptitud, n, rng, k=3, **_):
    """Torneo: se toman k cromosomas al azar y el de mejor puntaje es el padre."""
    aptitud = np.asarray(aptitud, float)
    k = max(1, int(k))
    competidores = rng.integers(0, aptitud.size, size=(n, k))
    gana = np.argmax(aptitud[competidores], axis=1)
    return competidores[np.arange(n), gana]


def seleccion_ranking(aptitud, n, rng, presion=1.7, **_):
    """Ranking lineal: la probabilidad depende del puesto en la tabla, no del puntaje.

        p(r) = [ (2 − s) + 2 (s − 1) r / (N − 1) ] / N,   r = 0 el peor, N − 1 el mejor
    """
    aptitud = np.asarray(aptitud, float)
    N = aptitud.size
    if N == 1:
        return np.zeros(n, dtype=int)
    s = float(np.clip(presion, 1.0, 2.0))
    rangos = np.empty(N)
    rangos[np.argsort(aptitud)] = np.arange(N)
    prob = ((2 - s) + 2 * (s - 1) * rangos / (N - 1)) / N
    return seleccion_ruleta(prob, n, rng)


SELECCIONES = {"Ruleta": seleccion_ruleta, "Torneo": seleccion_torneo, "Ranking": seleccion_ranking}


# ---------------------------------------------------------------------------
# CRUCE (cada función recibe dos padres y devuelve dos hijos)
# ---------------------------------------------------------------------------

def _dos_cortes(L, rng):
    a, b = np.sort(rng.choice(L + 1, 2, replace=False))
    if b - a < 1:
        b = a + 1
    return int(a), int(b)


def pmx(p1, p2, rng, cortes=None):
    """PMX, cruce parcialmente mapeado (Goldberg y Lingle, 1985).

    1. Se elige un tramo [a, b) y el hijo copia ese tramo del padre 1.
    2. Las posiciones de fuera se llenan con el padre 2; si el jugador ya
       está en el tramo, se sigue el "mapa" tramo_p1 -> tramo_p2 hasta dar
       con uno que no esté repetido.

        p1 = 3 1 | 4 2 | 5        tramo copiado de p1: 4 2   (mapa 4->1, 2->7…)
        p2 = 2 4 | 1 7 | 3
    """
    L = p1.size
    a, b = cortes if cortes is not None else _dos_cortes(L, rng)

    def hijo(x, y):
        h = np.full(L, -1)
        h[a:b] = x[a:b]
        en_tramo = set(x[a:b].tolist())
        mapa = {int(x[i]): int(y[i]) for i in range(a, b)}
        for i in list(range(0, a)) + list(range(b, L)):
            g = int(y[i])
            while g in en_tramo:
                g = mapa[g]
            h[i] = g
        return h
    return hijo(p1, p2), hijo(p2, p1)


def ox(p1, p2, rng, cortes=None):
    """OX, cruce de orden (Davis, 1985).

    El hijo copia el tramo [a, b) del padre 1 y el resto de jugadores se
    ponen en el ORDEN en que aparecen en el padre 2, empezando después del
    tramo y dando la vuelta. Conserva el orden relativo del otro padre.
    """
    L = p1.size
    a, b = cortes if cortes is not None else _dos_cortes(L, rng)

    def hijo(x, y):
        h = np.full(L, -1)
        h[a:b] = x[a:b]
        usados = set(x[a:b].tolist())
        orden = [int(y[(b + i) % L]) for i in range(L) if int(y[(b + i) % L]) not in usados]
        posiciones = [(b + i) % L for i in range(L - (b - a))]
        for pos, g in zip(posiciones, orden):
            h[pos] = g
        return h
    return hijo(p1, p2), hijo(p2, p1)


def ciclos(p1, p2, rng=None, **_):
    """CX, cruce por ciclos (Oliver, Smith y Holland, 1987).

    Se recorren los "ciclos" de posiciones: se arranca en la posición 0, se
    mira qué jugador tiene ahí el padre 2, se busca ese jugador en el padre 1
    y se salta a esa posición, hasta volver al inicio. Los ciclos se toman
    alternando de un padre y del otro. Cada jugador queda en una posición
    que ocupaba en alguno de sus padres: hereda POSICIONES, no orden.
    """
    L = p1.size
    h1 = np.full(L, -1)
    h2 = np.full(L, -1)
    pos_en_p1 = np.empty(L, int)
    pos_en_p1[p1] = np.arange(L)
    visitado = np.zeros(L, bool)
    numero_ciclo = 0
    for inicio in range(L):
        if visitado[inicio]:
            continue
        i = inicio
        while not visitado[i]:
            visitado[i] = True
            if numero_ciclo % 2 == 0:
                h1[i], h2[i] = p1[i], p2[i]
            else:
                h1[i], h2[i] = p2[i], p1[i]
            i = pos_en_p1[p2[i]]
        numero_ciclo += 1
    return h1, h2


CRUCES = {"PMX": pmx, "OX (orden)": ox, "Ciclos (CX)": ciclos}


def cruzar_poblacion(padres, pc, cruce, rng):
    """Cruza los padres de dos en dos (0-1, 2-3, …) con probabilidad pc."""
    hijos = padres.copy()
    funcion = CRUCES[cruce]
    for i in range(0, padres.shape[0] - 1, 2):
        if rng.random() < pc:
            hijos[i], hijos[i + 1] = funcion(padres[i], padres[i + 1], rng)
    return hijos


# ---------------------------------------------------------------------------
# MUTACIÓN (cada cromosoma muta con probabilidad pm)
# ---------------------------------------------------------------------------

def mut_intercambio(c, rng):
    """Intercambio: dos jugadores cambian de puesto (o un titular con uno del banco)."""
    c = c.copy()
    i, j = rng.choice(c.size, 2, replace=False)
    c[i], c[j] = c[j], c[i]
    return c


def mut_insercion(c, rng):
    """Inserción: un jugador sale de su lugar y entra en otro; los del medio se corren."""
    c = list(c)
    i, j = rng.choice(len(c), 2, replace=False)
    g = c.pop(i)
    c.insert(j, g)
    return np.array(c)


def mut_inversion(c, rng):
    """Inversión: un tramo del cromosoma se voltea al revés.

    El tramo tiene al menos dos genes: voltear un solo gen no cambiaría nada.
    """
    c = c.copy()
    a, b = np.sort(rng.choice(c.size + 1, 2, replace=False))
    while b - a < 2:
        a, b = np.sort(rng.choice(c.size + 1, 2, replace=False))
    c[a:b] = c[a:b][::-1]
    return c


MUTACIONES = {"Intercambio": mut_intercambio, "Inserción": mut_insercion, "Inversión": mut_inversion}


def mutar_poblacion(hijos, pm, mutacion, rng):
    funcion = MUTACIONES[mutacion]
    salida = hijos.copy()
    for i in range(hijos.shape[0]):
        if rng.random() < pm:
            salida[i] = funcion(hijos[i], rng)
    return salida


def es_permutacion(c, n):
    """True si c contiene cada número de 0 a n-1 exactamente una vez."""
    return np.array_equal(np.sort(c), np.arange(n))
