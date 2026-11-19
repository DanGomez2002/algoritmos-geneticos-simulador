"""
abrir_codigo.py
---------------
Clic derecho en el simulador -> se abre VS Code en la función que hace ese
trabajo. La línea se calcula en el momento del clic leyendo el archivo con
`ast`, así nunca queda desactualizada aunque se edite el código.
"""

import ast
import os
import shutil
import subprocess
import sys
import urllib.parse

CARPETA = os.path.dirname(os.path.abspath(__file__))


def ruta(archivo):
    """Ruta absoluta de un archivo del proyecto."""
    return os.path.join(CARPETA, archivo)


def linea_de(archivo, nombre, buscar=None):
    """Línea donde se define `nombre` en `archivo`.

    `nombre` puede ser "funcion" o "Clase.metodo". Si se da `buscar`, se
    devuelve la primera línea dentro de la función que contiene ese texto.
    Si no se encuentra, devuelve 1.
    """
    with open(ruta(archivo), encoding="utf-8") as fh:
        texto = fh.read()
    arbol = ast.parse(texto)
    partes = nombre.split(".")

    def buscar_en(nodos, resto):
        for nodo in nodos:
            if isinstance(nodo, (ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef)) \
                    and nodo.name == resto[0]:
                if len(resto) == 1:
                    return nodo
                return buscar_en(nodo.body, resto[1:])
        return None

    nodo = buscar_en(arbol.body, partes)
    if nodo is None:
        # También puede ser una variable global (por ejemplo, INGREDIENTES).
        for n in arbol.body:
            if isinstance(n, ast.Assign) and any(
                    isinstance(t, ast.Name) and t.id == nombre for t in n.targets):
                return n.lineno
        return 1
    if buscar:
        lineas = texto.splitlines()
        fin = getattr(nodo, "end_lineno", len(lineas))
        for i in range(nodo.lineno - 1, fin):
            if buscar in lineas[i]:
                return i + 1
    return nodo.lineno


def abrir(archivo, nombre, buscar=None):
    """Abre VS Code en archivo:línea sin bloquear el simulador.

    Devuelve (abierto, "archivo:línea"). Nunca lanza excepción.
    """
    try:
        linea = linea_de(archivo, nombre, buscar)
    except Exception:
        linea = 1
    destino = f"{archivo}:{linea}"
    camino = ruta(archivo)
    try:
        code = shutil.which("code") or shutil.which("code.cmd")
        if code:
            opciones = {}
            if sys.platform.startswith("win"):
                opciones["creationflags"] = getattr(subprocess, "CREATE_NO_WINDOW", 0)
            subprocess.Popen([code, "-g", f"{camino}:{linea}"], **opciones)
            return True, destino
        url = "vscode://file/" + urllib.parse.quote(camino.replace("\\", "/"), safe="/:") + f":{linea}"
        if sys.platform.startswith("win"):
            os.startfile(url)                   # noqa: S606 (solo Windows)
            return True, destino
        comando = "open" if sys.platform == "darwin" else "xdg-open"
        if shutil.which(comando):
            subprocess.Popen([comando, url])
            return True, destino
    except Exception:
        pass
    return False, destino
