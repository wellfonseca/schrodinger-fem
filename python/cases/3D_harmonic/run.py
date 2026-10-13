#!/usr/bin/env python3
"""Case 3D — oscilador harmonic isotropico.  Analytic: E_n = 2.sqrt(a)(n+3/2)"""
import sys, os
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', '..'))   # python/
sys.path.insert(0, os.path.join(_here, '..'))          # python/cases/
import numpy as np
from _common import run_case

L, n, alpha = 12.0, 14, 0.25
# Levels: E = 2.sqrt(alpha).(N + d/2), N = 0,1,2,...  with degeneracy (N+1)(N+2)/2
_sd = 2.0 * np.sqrt(alpha)
LEVELS = []
for _N in range(0, 6):
    LEVELS += [_sd * (_N + 1.5)] * ((_N + 1) * (_N + 2) // 2)   # degeneracy em 3D

run_case(3, L, n, alpha, "3D harmonic V=0.25(x^2+y^2+z^2) (L=12, 14x14x14)",
     ana_fn=lambda k: LEVELS[k] if k < len(LEVELS) else None)
