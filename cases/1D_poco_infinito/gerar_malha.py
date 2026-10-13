#!/usr/bin/env python3
"""Gera a malha Elmer do caso 1D_poco_infinito (1, L=40.0, n=250)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts'))
from gera_malha_elmer import gera_1d, gera_2d, gera_3d

d = os.path.dirname(os.path.abspath(__file__))
os.chdir(d)
gera_1d('malha', 40.0, 250, -20.0)
