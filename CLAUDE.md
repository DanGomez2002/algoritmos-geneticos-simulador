# Contexto del proyecto (para Claude Code)

Entrega final de Computacional I (profesor Jesús Emilio Pinto Lopera), autor Carlos Daniel Gómez Murcia. Tema: algoritmos genéticos con codificación por permutación. Caso: alineación de un equipo de fútbol (18 jugadores, 11 puestos + banco). Responde y comenta el código en español.

## Lo que pide el profesor (clase del 18-19 nov 2026)
* Ejercicio de permutación (no rutas ni repartos), aprobado por correo antes de entregarlo.
* Simulación que muestre cómo se organizan los genes generación tras generación (la mejor/élite) y cómo queda al final. Unity recomendado, no obligatorio.
* Gráfica del mejor y del promedio por generación.
* Experimentos: tamaño de población (10, 50, 200…), generaciones necesarias, varias corridas.
* En la sustentación preguntará por la población, los cromosomas, los genes, la función de evaluación, el elitismo, el cruce, la mutación y las condiciones de parada.
* Sustentaciones desde el lunes 30 de noviembre.

## Estado (19 nov 2026)
* Simulador terminado; `pruebas.py` tiene 21 comprobaciones y todas pasan.
* Resultados en `resultados/` (`generar_resultados.py`, unos seis minutos). Incluyen experimentos con el algoritmo y con el fenómeno: química contra método húngaro, barrido del peso de la química (λ = 0, 1, 3, 5, 10), formaciones e importancia de cada titular (se lesiona uno por uno).
* La tabla de genes dibuja cada número como una figura vectorial reutilizada (18 objetos en vez de 540 textos); no volver a `ax.text` por celda, porque la ventana se vuelve lenta.
* La mutación por inversión exige un tramo de al menos 2 genes; con eso PMX + inversión da 30 % de éxito y ciclos + inversión 10 %.
* Un extremo (ED/EI) jugando de volante por banda no se marca como fuera de posición.
* Algunas alineaciones empatan en puntaje (dos centrales o dos volantes que cambian de lado, Ospina y Calderón con 71,4 como volante central); entre equipos puede salir una u otra con el mismo puntaje.
* El profesor aprobó el ejercicio «más o menos»; pidió terminar de definir las características del fenómeno de estudio.
* Pendiente: informe Word, app de estudio y guion de sustentación.

## Reglas
* Al cambiar algo, correr `python pruebas.py` y `python generar_resultados.py`, y avisar qué cifras del informe cambian.
* No mencionar git en la guía de ejecución: se entrega en .zip.
