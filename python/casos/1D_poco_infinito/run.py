#!/usr/bin/env python3
"""Caso 1D — poco infinito.  Analitico: E_n = (n.pi/L)^2"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from _comum import roda

L, n = 40.0, 250
roda(1, L, n, 0.0, "1D poco infinito (L=40, 250 elementos)",
     ana_fn=lambda k: ((k + 1) * np.pi / L)**2)
