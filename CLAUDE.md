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
* Simulador terminado y 19 comprobaciones en `pruebas.py` que pasan (las 2 nuevas: portero lesionado no cuenta como suplente e inversión siempre cambia algo; además la semilla negativa ya no falla en la ventana).
* Resultados en `resultados/` (`generar_resultados.py`).
* Pendiente: aprobación del profesor, informe Word, app de estudio y guion de sustentación.

## Reglas
* Al cambiar algo, correr `python pruebas.py` y `python generar_resultados.py`, y avisar qué cifras del informe cambian.
* No mencionar git en la guía de ejecución: se entrega en .zip.
