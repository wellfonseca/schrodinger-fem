#!/usr/bin/env python3
"""Gera a malha Elmer do caso 3D_harmonico (3, L=12.0, n=14)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts'))
from gera_malha_elmer import gera_1d, gera_2d, gera_3d

d = os.path.dirname(os.path.abspath(__file__))
os.chdir(d)
gera_3d('malha', 12.0, 12.0, 12.0, 14, 14, 14, -6.0, -6.0, -6.0)
