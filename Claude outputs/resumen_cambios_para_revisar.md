# Cambios al simulador del algoritmo genético: ¿son válidos?

## Contexto
Proyecto académico en Python (numpy, matplotlib, scipy). Un algoritmo genético con codificación por permutación arma la alineación de un equipo de fútbol:
* **Cromosoma:** permutación de 18 jugadores. Los genes 0 a 10 son los 11 puestos y los genes 11 a 17 son el banco.
* **Puntaje:** Σ rendimiento(jugador, puesto) + λ·química entre puestos vecinos − 40 si no hay portero suplente en el banco − 300 por cada lesionado titular.
* **Operadores:** selección por ruleta, torneo o ranking; cruces PMX, OX o CX; mutación por intercambio, inserción o inversión; elitismo.
* **Ventana:** interfaz interactiva hecha con widgets de matplotlib.

Pedí revisar todo el proyecto. Primero se escribió una batería de 23 verificaciones independientes:
* el puntaje recalculado con una fórmula escrita a mano;
* propiedades de PMX, OX y CX con todos los puntos de corte;
* frecuencias de selección contra las fórmulas teóricas;
* elitismo y condiciones de parada exactas;
* las 81 combinaciones de operadores;
* óptimo exacto del método húngaro cuando no hay química;
* búsqueda local independiente que no logra superar al algoritmo genético;
* entradas raras en la interfaz.

Todo pasó, pero salieron 3 detalles. Después se analizó la lentitud de la ventana.

## Cambio 1. Semilla negativa en la ventana (interfaz.py)
**Problema:** escribir `-3` en la caja "Semilla" lanzaba `ValueError` en `np.random.default_rng(-3)` y la corrida no se reiniciaba.

**Cambio:**
```python
semilla=abs(int(self._numero(self.t_sem, 7)))   # NumPy no acepta semillas negativas
```

## Cambio 2. Un portero lesionado no cuenta como suplente (evaluacion.py)
**Problema:** la penalización de 40 por "banco sin portero suplente" se evitaba con un portero **lesionado** en el banco.

**Antes:**
```python
sin_portero = ~pl.ES_PORTERO[banco].any(axis=1)
```

**Ahora:**
```python
portero_sano = pl.ES_PORTERO[banco] & ~np.isin(banco, self.lesionados)
sin_portero = ~portero_sano.any(axis=1)
```

El método húngaro de referencia (`optimo_sin_quimica`) no modela esta penalización. Antes tampoco lo hacía: solo es una referencia sin química.

## Cambio 3. La mutación por inversión siempre cambia algo (operadores.py)
**Problema:** se elegían dos cortes `a < b` en `0..L` sin repetir. Cuando `b − a = 1`, se "invertía" un solo gen y no cambiaba nada. Pasaba en 18 de los 171 pares de cortes, cerca del 11 % de las mutaciones.

**Ahora:** se vuelve a sortear hasta que el tramo tenga al menos 2 genes, lo que mantiene uniforme la distribución entre los tramos válidos.
```python
a, b = np.sort(rng.choice(c.size + 1, 2, replace=False))
while b - a < 2:
    a, b = np.sort(rng.choice(c.size + 1, 2, replace=False))
c[a:b] = c[a:b][::-1]
```

**Efecto en los resultados (30 semillas, población 50, 300 generaciones):** solo cambian las 3 filas de la tabla de operadores que usan inversión.

| Combinación | Puntaje medio | Desviación | Éxito |
|---|---|---|---|
| PMX + Inversión | 844,6 → 846,6 | 6,18 → 4,73 | 23,3 % → 30,0 % |
| OX + Inversión | 835,4 → 835,0 | 8,23 → 8,22 | sin cambio |
| CX + Inversión | 839,0 → 839,9 | 6,41 → 6,34 | 3,3 % → 10,0 % |

## Cambio 4. Velocidad de la ventana (graficas.py)
**Problema:** el algoritmo tarda de 1 a 10 ms por generación, pero redibujar la ventana tardaba 760 ms por generación en la ventana real (TkAgg).

**Diagnóstico:** se midió el tiempo de dibujo de cada panel. El panel "genes de la población" se llevaba cerca de la mitad, porque ponía un número en cada celda con `ax.text`: 30 filas × 18 genes = 540 objetos Text. Matplotlib tarda cerca de 1 ms en dibujar cada texto.

**Cambio:** cada número del 1 al 18 se convierte una sola vez en un `TextPath` centrado y se estampa en todas sus celdas con un `PathCollection`. Así son 18 objetos en vez de 540.
```python
escala = Affine2D().scale(1 / 72) + ax.figure.dpi_scale_trans   # puntos -> píxeles, respeta el dpi al guardar
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
```

**Mediciones:**

| Medición | Antes | Ahora |
|---|---|---|
| Dibujo del panel de genes | 489 ms | 45 ms |
| Redibujo completo (backend Agg) | 1373 ms | 540 ms |
| Ventana real (TkAgg) | 759 ms por generación (1,3 gen/s) | ~500 ms por generación (2 gen/s) |

**Comprobación visual:** las imágenes del panel antes y después se compararon a 100 y a 150 ppp, y se ven iguales.

**Alternativas descartadas:**
* **Dibujar solo cada N generaciones:** acelera, pero esconde justo lo que el profesor quiere ver, que es cómo se organizan los genes generación tras generación.
* **Blitting o reutilizar objetos de matplotlib:** el tiempo que queda está repartido entre los paneles (de 50 a 90 ms cada uno), así que daría poca mejora a cambio de más complejidad. El estudiante tiene que poder explicar el código en la sustentación.

## Cómo se verificó
* `python pruebas.py`: pasa 19 de 19. Se agregaron 2 pruebas: el portero lesionado no cuenta como suplente y la inversión siempre cambia algo. La prueba de la interfaz ahora también cubre la semilla negativa.
* Las 23 verificaciones independientes pasan. La fórmula escrita a mano se actualizó a la regla nueva del portero sano.
* `generar_resultados.py` es reproducible: dos corridas con el mismo código dan un JSON idéntico byte a byte. El cambio 4 no alteró ningún número.

## Preguntas para la revisión
1. ¿Es correcto conceptualmente que un portero lesionado no cuente como suplente, o habría que modelarlo de otra forma?
2. ¿Sigue siendo una "mutación por inversión" válida si se exige un tramo de al menos 2 genes? ¿Es preferible el muestreo por rechazo o elegir directamente entre los pares válidos?
3. ¿`abs()` es la mejor forma de tratar una semilla negativa, o sería mejor rechazarla y avisar al usuario?
4. ¿El uso de `TextPath` y `PathCollection` con `dpi_scale_trans` es correcto y robusto, por ejemplo al cambiar el tamaño de la ventana o en pantallas HiDPI?
5. ¿Se pasa por alto alguna otra optimización simple de la ventana?
