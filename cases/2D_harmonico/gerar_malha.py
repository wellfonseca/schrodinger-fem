#!/usr/bin/env python3
"""Gera a malha Elmer do caso 2D_harmonico (2, L=20.0, n=40)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts'))
from gera_malha_elmer import gera_1d, gera_2d, gera_3d

d = os.path.dirname(os.path.abspath(__file__))
os.chdir(d)
gera_2d('malha', 20.0, 20.0, 40, 40, -10.0, -10.0)
