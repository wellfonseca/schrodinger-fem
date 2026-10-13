#!/usr/bin/env python3
"""Case 1D — infinite well.  Analytic: E_n = (n.pi/L)^2"""
import sys, os
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, '..', '..'))   # python/
sys.path.insert(0, os.path.join(_here, '..'))          # python/cases/
import numpy as np
from _common import run_case

L, n = 40.0, 250
run_case(1, L, n, 0.0, "1D infinite well (L=40, 250 elements)",
     ana_fn=lambda k: ((k + 1) * np.pi / L)**2)
