#!/usr/bin/env python3
"""Generate the Elmer mesh for case 1D_harmonic (1, L=20.0, n=400)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from generate_elmer_mesh import generate_1d, generate_2d, generate_3d

d = os.path.dirname(os.path.abspath(__file__))
os.chdir(d)
generate_1d('mesh', 20.0, 400, -10.0)
