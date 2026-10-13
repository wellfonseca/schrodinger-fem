#!/usr/bin/env python3
"""Gera a malha Elmer do caso 1D_morse (1, L=10.0, n=400)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts'))
from gera_malha_elmer import gera_1d, gera_2d, gera_3d

d = os.path.dirname(os.path.abspath(__file__))
os.chdir(d)
gera_1d('malha', 10.0, 400, 0.0)
