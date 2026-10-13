#!/usr/bin/env python3
"""Case 3D — cubic well.  Analytic: E = pi^2(nx^2+ny^2+nz^2)/L^2"""
import sys, os
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', '..'))   # python/
sys.path.insert(0, os.path.join(_here, '..'))          # python/cases/
import numpy as np
from _common import run_case

L, n = 6.0, 10
def ana(k):
    tri = [(1,1,1), (1,1,2), (1,2,1), (2,1,1), (2,2,1)]  # E ~ nx^2+ny^2+nz^2
    if k >= len(tri): return None
    a, b, c = tri[k]
    return np.pi**2 * (a*a + b*b + c*c) / L**2

run_case(3, L, n, 0.0, "3D cubic well (L=6, 10x10x10)", ana)
