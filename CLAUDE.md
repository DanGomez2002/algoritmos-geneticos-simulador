# Contexto del proyecto (para Claude Code)

Entrega final de Computacional I (profesor Jesús Emilio Pinto Lopera), autor Carlos Daniel Gómez Murcia. Tema: algoritmos genéticos con codificación por permutación. Caso: alineación de un equipo de fútbol (18 jugadores, 11 puestos + banco). Responde y comenta el código en español.

## Lo que pide el profesor (clase del 18-19 nov 2026)
* Ejercicio de permutación (no rutas ni repartos), aprobado por correo antes de entregarlo.
* Simulación que muestre cómo se organizan los genes generación tras generación (la mejor/élite) y cómo queda al final. Unity recomendado, no obligatorio.
* Gráfica del mejor y del promedio por generación.
* Experimentos: tamaño de población (10, 50, 200…), generaciones necesarias, varias corridas.
* En la sustentación preguntará por la población, los cromosomas, los genes, la función de evaluación, el elitismo, el cruce, la mutación y las condiciones de parada.
* Sustentaciones desde el lunes 30 de noviembre.

## Estado (26 nov 2026)
* Simulador terminado; `pruebas.py` tiene 23 comprobaciones y todas pasan.
* 25 nov: la ventana trae activos los tres criterios de parada: máximo 200 generaciones, 50 sin mejora y puntaje objetivo = mejor conocido de la formación (`OBJETIVO_CONOCIDO` en `interfaz.py`; se apaga con λ ≠ 3 o con lesionados). La población por defecto de la ventana es 100. Así suele parar entre la generación 25 y 155 en vez de llegar a 200. `Configuracion` y `experimentos.py` no cambiaron, así que ninguna cifra de `resultados/` cambia.
* 25 nov: la tabla de genes muestra los 15 mejores («Genes de la población, generación N (15 mejores de M)», etiquetas «1.º» y «15.º», leyenda debajo) y `ajustar_a_pantalla` en `interfaz.py` reduce el dpi de la figura cuando la pantalla es pequeña (probado en 1366 × 768 y 1920 × 1080 con TkAgg). 26 nov: en pantalla pequeña ya no encoge todo igual (las letras quedaban de 6-7 px en el portátil de Daniel, 1366 × 768): maximiza la ventana, quita la barra de herramientas y aplica una lupa `AUMENTO_LETRA = 1.25` (más de eso choca). `_ajustar_letra` mide el panel de texto y los titulares y baja la letra hasta que quepan.
* Resultados en `resultados/` (`generar_resultados.py`, unos siete minutos). Incluyen experimentos con el algoritmo (población, operadores, elitismo, mutación, selección, criterios de parada) y con el fenómeno: química contra método húngaro, barrido del peso de la química (λ = 0, 1, 3, 5, 10), formaciones e importancia de cada titular (se lesiona uno por uno).
* 26 nov: las cifras cambiaban entre Windows y Linux porque alineaciones empatadas (mismo puntaje) se desempataban distinto por diferencias de redondeo en el último decimal. Se corrigió en `algoritmo.py` (puntaje redondeado a 6 decimales y orden estable) y en `operadores.py` (ranking con orden estable). Con eso los resultados deben ser idénticos en cualquier computador; se regeneró todo y se agregó el experimento de criterios de parada (`criterios_parada.png`). Caso base (pob. 50, ruleta, PMX + intercambio): éxito 83,3 %, promedio 850,56, desv. 1,66, peor 844,46, llegada media 127,9 (máx. 299). Estas son las cifras del informe.
* La tabla de genes dibuja cada número como una figura vectorial reutilizada (18 objetos en vez de 540 textos); no volver a `ax.text` por celda, porque la ventana se vuelve lenta.
* La mutación por inversión exige un tramo de al menos 2 genes.
* Un extremo (ED/EI) jugando de volante por banda no se marca como fuera de posición.
* Algunas alineaciones empatan en puntaje (dos centrales o dos volantes que cambian de lado, Ospina y Calderón con 71,4 como volante central); entre equipos puede salir una u otra con el mismo puntaje.
* El profesor aprobó el ejercicio «más o menos»; pidió terminar de definir las características del fenómeno de estudio.
* Informe Word: primera versión completa (26 nov), 65 páginas con la plantilla de la entrega 3 (`Informe_Algoritmos_Geneticos_Carlos_Gomez_Alvaro_Garcia.docx`), autores Carlos Daniel Gómez Murcia y Alvaro Javier Garcia Montaña. 16 tablas, 14 figuras (incluidas las ventanas «Genes» y «Resumen»), 5 listados y ecuaciones F1 a F10 conectadas con las funciones. Las cifras salen de `resultados/resultados.json`: falta confirmar que la corrida en Windows con el código del 26 nov da el mismo JSON (la de las 8:34 era del código viejo). Corrida en Windows de las 9:02 a las 9:10 (`algoritmo.py` y `operadores.py` del 26 nov, 8:53): el caso base da exactamente 83,3 %, 850,56, 1,66, 844,46, 127,9 y 299; no se comparó el JSON completo contra el de Linux.
* 26 nov: `figura_interfaz` captura la ventana cuando se cumple el criterio de parada (gen. 42), las etiquetas de `criterios_parada.png` van en dos líneas y la tabla de «Ver plantilla» ya no corta los nombres. En el PC de Daniel se movió la lista de titulares debajo de «Mutación» y se compactó el panel de texto.
* 26 nov: la corrida en Windows de las 9:10 da un `resultados.json` idéntico al de Linux (cero diferencias): las cifras del informe quedan confirmadas.
* 26 nov: el botón «Guardar imagen» se cambió por tres botones en la misma fila: «Plantilla», «Genes» y «Resumen». «Genes» (`Simulador._ver_genes`, `graficas.dibujar_mejores`) muestra el mejor cromosoma de cada generación con mejora, con un deslizador «Hasta la generación». «Resumen» (`Simulador._ver_resumen`, `graficas.fig_resumen`, `cambios_entre`) muestra evolución con banda peor-mejor, diversidad, tabla de mejoras con lo que cambió en los titulares (sin contar cambios entre puestos iguales) y datos de la corrida. Ambas tienen «Guardar imagen» (`_ventana_con_guardar`), que guarda la principal y la ventana abierta. Se probó y se descartó mostrar toda la población en la ventana de genes. «Plantilla» pone los pesos de los roles debajo de los jugadores; las tres ventanas usan `ajustar_a_pantalla`. Prueba `test_ventanas_genes_y_resumen` (23 en total). `generar_resultados.py` guarda `ventana_genes.png` y `resumen.png`.
* Pendiente: app de estudio y guion de sustentación.

## Reglas
* Al cambiar algo, correr `python pruebas.py` y `python generar_resultados.py`, y avisar qué cifras del informe cambian.
* No mencionar git en la guía de ejecución: se entrega en .zip.
