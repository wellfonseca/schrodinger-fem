#!/usr/bin/env python3
"""Caso 2D — poco quadrado.  Analitico: E = pi^2(nx^2+ny^2)/L^2"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from _comum import roda

L, n = 10.0, 40
def ana(k):
    # niveis ordenados: 1,2,2,4,5  ->  (nx,ny) = (1,1);(1,2);(2,1);(2,2);(1,3)
    pares = [(1,1), (1,2), (2,1), (2,2), (1,3)]
    if k >= len(pares): return None
    a, b = pares[k]
    return np.pi**2 * (a*a + b*b) / L**2

roda(2, L, n, 0.0, "2D poco quadrado (L=10, 40x40)", ana)
