"""
main.py
-------
    python main.py                        abre el simulador
    python main.py --captura prueba.png   corre 40 generaciones sin ventana y guarda la imagen
"""

import argparse
import sys


def principal():
    parser = argparse.ArgumentParser(description="Alineación de fútbol con algoritmo genético")
    parser.add_argument("--captura", metavar="ARCHIVO.png")
    parser.add_argument("--generaciones", type=int, default=40)
    args = parser.parse_args()
    if args.captura:
        import matplotlib
        matplotlib.use("Agg")
        from interfaz import Simulador
        sim = Simulador()
        for _ in range(args.generaciones):
            if sim.ag.terminado():
                break
            sim.ag.paso()
        sim._dibujar()
        sim.fig.savefig(args.captura, dpi=100)
        print(f"Captura guardada en {args.captura} (generación {sim.ag.generacion})")
        return 0
    from interfaz import lanzar
    lanzar()
    return 0


if __name__ == "__main__":
    sys.exit(principal())
