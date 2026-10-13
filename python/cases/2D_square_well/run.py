#!/usr/bin/env python3
"""Case 2D — square well.  Analytic: E = pi^2(nx^2+ny^2)/L^2"""
import sys, os
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', '..'))   # python/
sys.path.insert(0, os.path.join(_here, '..'))          # python/cases/
import numpy as np
from _common import run_case

L, n = 10.0, 40
def ana(k):
    # niveis ordenados: 1,2,2,4,5  ->  (nx,ny) = (1,1);(1,2);(2,1);(2,2);(1,3)
    pares = [(1,1), (1,2), (2,1), (2,2), (1,3)]
    if k >= len(pares): return None
    a, b = pares[k]
    return np.pi**2 * (a*a + b*b) / L**2

run_case(2, L, n, 0.0, "2D square well (L=10, 40x40)", ana)
