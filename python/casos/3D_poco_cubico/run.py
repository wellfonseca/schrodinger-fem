#!/usr/bin/env python3
"""Caso 3D — poco cubico.  Analitico: E = pi^2(nx^2+ny^2+nz^2)/L^2"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from _comum import roda

L, n = 6.0, 10
def ana(k):
    tri = [(1,1,1), (1,1,2), (1,2,1), (2,1,1), (2,2,1)]  # E ~ nx^2+ny^2+nz^2
    if k >= len(tri): return None
    a, b, c = tri[k]
    return np.pi**2 * (a*a + b*b + c*c) / L**2

roda(3, L, n, 0.0, "3D poco cubico (L=6, 10x10x10)", ana)
