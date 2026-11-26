# Alineación de un equipo de fútbol con algoritmo genético

Entrega final de Computacional I. Autor: Carlos Daniel Gómez Murcia, Universidad de la Amazonia.

Un entrenador tiene una plantilla de 18 jugadores y debe decidir quién juega en cada uno de los 11 puestos de su formación y quiénes van al banco. El algoritmo genético busca la alineación con el mayor puntaje, teniendo en cuenta lo siguiente:
* el rendimiento de cada jugador en cada puesto, calculado con sus estadísticas y el perfil del rol;
* la química entre jugadores que juegan cerca;
* el pie hábil en las bandas;
* el estado físico;
* que haya un portero suplente;
* que no jueguen los lesionados.

## Cómo ejecutarlo

1. Tener Python 3.10 o superior.
2. En una terminal abierta en esta carpeta: `pip install -r requirements.txt`
3. Abrir el simulador: `python main.py`
4. Comprobar que todo funciona: `python pruebas.py` (veintitrés comprobaciones).
5. Regenerar las figuras y los números del informe: `python generar_resultados.py` (unos siete minutos; quedan en `resultados/`).

## Cómo está planteado

| Concepto | En este ejercicio |
|---|---|
| Cromosoma | Una alineación: permutación de los 18 jugadores. Los genes 1 a 11 son los puestos de la formación y los genes 12 a 18 son el banco |
| Gen | El jugador que ocupa una posición del cromosoma; trae sus estadísticas (ritmo, tiro, pase, regate, defensa, físico, portería), su pie y su estado físico |
| Población | Conjunto de alineaciones (100 por defecto en la ventana; los experimentos usan 50) |
| Evaluación | Σ rendimiento de los titulares + λ · química entre vecinos − penalizaciones |
| Selección | Ruleta (por defecto), torneo o ranking |
| Cruce | PMX, OX (orden) o ciclos (CX), que nunca repiten jugadores |
| Mutación | Intercambio (por defecto), inserción o inversión |
| Elitismo | Porcentaje de las mejores alineaciones que pasa intacto (10 % por defecto) |
| Parada | El primero que se cumpla: máximo de generaciones (200), generaciones sin mejora (50) o puntaje objetivo (el mejor conocido de la formación: 851,2 en el 4-4-2, 848,2 en el 4-3-3 y 845,5 en el 3-5-2; se apaga con otro λ o con lesionados) |

## Archivos

| Archivo | Qué contiene |
|---|---|
| `plantilla.py` | Jugadores y estadísticas, pesos de cada rol, formaciones y química |
| `evaluacion.py` | Matriz de rendimiento, función de evaluación, penalizaciones y referencia exacta sin química (método húngaro) |
| `operadores.py` | Selección, cruces PMX, OX y ciclos, mutaciones |
| `algoritmo.py` | El ciclo del algoritmo genético, el elitismo y las condiciones de parada |
| `graficas.py` | Cancha, evolución del puntaje, genes de la población, matriz de rendimiento y tabla de la plantilla |
| `interfaz.py` | La ventana del simulador |
| `experimentos.py` | Experimentos con el algoritmo (población, operadores, elitismo, mutación, selección; 30 semillas) y con el fenómeno (química, formaciones, importancia de cada titular al lesionarlo) |
| `generar_resultados.py` | Figuras y `resultados/resultados.json` |
| `abrir_codigo.py` | Clic derecho que abre el código en VS Code |
| `pruebas.py` | Comprobaciones automáticas |

## Uso de la ventana

* **Ejecutar / Pausa** corre una generación tras otra y **Un paso** avanza una sola. **Reiniciar** crea una población nueva con la semilla escrita.
* Selección, cruce, mutación, probabilidades, elitismo, tamaño del torneo y condiciones de parada se cambian **mientras corre**, sin perder la búsqueda.
* Cambiar la formación, la población, el peso de la química, la semilla o los lesionados empieza de nuevo.
* **Lesionados (#)** recibe los números de los jugadores separados por coma (por ejemplo `10, 3`).
* Si la pantalla es más pequeña que la ventana (por ejemplo, un portátil de 1366 × 768), la ventana se abre maximizada y sin la barra de herramientas de Matplotlib, y las letras se agrandan un 25 % para que se lean bien. Los textos de resultados ajustan su tamaño solos para caber completos.
* La tabla de genes muestra los 15 mejores cromosomas de la población, con el puesto en el ranking y el puntaje del primero y del último.
* **Ver generación** permite volver a cualquier generación para ver la mejor alineación y los genes de la población en ese momento.
* **Repetir evolución** reproduce toda la búsqueda desde la generación 0.
* **Plantilla** muestra las estadísticas de los jugadores y, debajo, los pesos de cada rol.
* **Genes** abre una ventana con el mejor cromosoma de cada generación en que mejoró la alineación (una fila por mejora, con su generación y su puntaje); el deslizador elige hasta qué generación se muestran.
* **Resumen** abre una ventana con la evolución del puntaje (mejor, promedio y rango de la población), la diversidad, la tabla de las mejoras con lo que cambió en los titulares y los datos de la corrida.
* En **Genes** y **Resumen**, el botón **Guardar imagen** guarda la ventana principal y la que está abierta (`captura_generacion_N.png` y `genes_generacion_N.png` o `resumen_generacion_N.png`).
* En la cancha, el anillo naranja punteado marca a un jugador que juega fuera de su línea natural. Las líneas verdes y rojas muestran la química entre vecinos.

## Clic derecho: ir al código

Clic derecho sobre cualquier control o gráfica abre VS Code en la función que hace ese trabajo. **Mayús + clic derecho** abre la función secundaria. En las listas de selección, cruce y mutación, el clic sobre una opción abre la función de esa opción. Para que funcione, el comando `code` debe estar en el PATH (en VS Code: *Shell Command: Install 'code' command in PATH*) o VS Code debe estar asociado a los enlaces `vscode://`.
