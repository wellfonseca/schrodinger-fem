#!/usr/bin/env python3
"""Generate the Elmer mesh for case 3D_harmonic (3, L=12.0, n=14)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from generate_elmer_mesh import generate_1d, generate_2d, generate_3d

d = os.path.dirname(os.path.abspath(__file__))
os.chdir(d)
generate_3d('mesh', 12.0, 12.0, 12.0, 14, 14, 14, -6.0, -6.0, -6.0)
