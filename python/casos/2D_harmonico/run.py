#!/usr/bin/env python3
"""Caso 2D — oscilador harmonico isotropico.  Analitico: E_n = 2.sqrt(a)(n+d/2)"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from _comum import roda

L, n, alpha = 30.0, 60, 0.25
# Niveis: E = 2.sqrt(alpha).(N + d/2), N = 0,1,2,...  com degenerescencia (N+1)
_sd = 2.0 * np.sqrt(alpha)
NIVEIS = []
for _N in range(0, 6):
    NIVEIS += [_sd * (_N + 1.0)] * (_N + 1)   # degenerescencia N+1 em 2D

roda(2, L, n, alpha, "2D harmonico V=0.25(x^2+y^2) (L=30, 60x60)",
     ana_fn=lambda k: NIVEIS[k] if k < len(NIVEIS) else None)
