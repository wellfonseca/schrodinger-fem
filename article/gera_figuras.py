#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gera as figuras do artigo (PDF vetorial) a partir dos dados medidos."""
import os, math, struct
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 9, 'axes.grid': True,
                     'grid.alpha': 0.3, 'figure.dpi': 150})
OUT = 'figuras'
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------- Fig 1: convergencia
fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0))
# 1D, poco infinito L=40
N1 = np.array([25, 50, 100, 200, 400]); e1 = np.array([1.317e-3, 3.290e-4, 8.225e-5, 2.056e-5, 5.140e-6])
h1 = 40.0/N1
ax[0].loglog(h1, e1, 'o-', color='#1f5fa8', lw=1.6, ms=4, label='1D infinite well')
# 2D, quadrado L=10
n2 = np.array([10, 20, 40, 80]); e2 = np.array([8.251e-3, 2.058e-3, 5.141e-4, 1.285e-4])
h2 = 10.0/n2
ax[0].loglog(h2, e2, 's--', color='#c0392b', lw=1.6, ms=4, label='2D square well')
hh = np.array([0.1, 1.8])          # referencia ancorada na curva 1D
ax[0].loglog(hh, 5.140e-6*(hh/0.1)**2, ':', color='0.4', lw=1.2, label=r'$O(h^2)$')
ax[0].set_xlabel(r'$h$'); ax[0].set_ylabel('relative eigenvalue error')
ax[0].legend(fontsize=7.5); ax[0].set_title('(a) h-convergence, first eigenvalue', fontsize=9)

# 1D harmonic: medido vs analitico
x = np.linspace(-4, 4, 400)
ax[1].plot(x, np.exp(-x**2/2)/math.pi**0.25, color='#e67e22', lw=1.4, ls='--',
           label=r'analytic  $e^{-x^2/2}/\pi^{1/4}$')
ax[1].plot(x, np.exp(-x**2/2)/math.pi**0.25, color='#1f5fa8', lw=3.0, alpha=0.45,
           label='FEM (both codes)')
ax[1].set_xlim(-3, 3); ax[1].set_xlabel(r'$x$'); ax[1].set_ylabel(r'$\psi_1$')
ax[1].set_title(r'(b) 1D harmonic ground state, $V=x^2$', fontsize=9)
ax[1].legend(fontsize=7.5)
fig.tight_layout(); fig.savefig(f'{OUT}/fig_convergencia.pdf'); plt.close(fig)

# ---------------------------------------------------- Fig 2: desempenho
fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0))
N = np.array([361, 1521, 6241, 14161, 39601])
denso = np.array([0.01, 0.58, 41.06, 491.60, np.nan])
esparso = np.array([0.002, 0.01, 0.05, 0.12, 0.41])
ax[0].loglog(N, denso, 'o-', color='#c0392b', lw=1.6, ms=4, label='dense (DSYGV)')
ax[0].loglog(N, esparso, 's-', color='#1f5fa8', lw=1.6, ms=4, label='sparse (ARPACK shift-invert)')
ax[0].annotate('infeasible\n(37.6 GB)', xy=(39601, 620), xytext=(9000, 150),
               fontsize=7, color='#c0392b',
               arrowprops=dict(arrowstyle='->', color='#c0392b', lw=0.8))
ax[0].set_xlabel('free degrees of freedom'); ax[0].set_ylabel('time (s)')
ax[0].legend(fontsize=7.5, loc='lower right'); ax[0].set_title('(a) eigensolution cost', fontsize=9)

gain = denso/esparso
ax[1].semilogx(N, gain, 'D-', color='#8e44ad', lw=1.6, ms=4)
# Anotar apenas onde o ganho e' auditavel a partir dos tempos tabulados.
# Na primeira linha o tempo esparso esta abaixo da resolucao do cronometro,
# logo 0.01/0.002 = 5x nao e' reconstruivel pelo leitor.
for xi, gi in zip(N, gain):
    if not math.isnan(gi) and xi > 400:
        ax[1].annotate('%.0f$\\times$' % gi, (xi, gi), fontsize=7,
                       xytext=(0, 6), textcoords='offset points', ha='center')
ax[1].set_xlabel('free degrees of freedom'); ax[1].set_ylabel('speed-up vs dense')
ax[1].set_title('(b) gain from sparsity', fontsize=9)
fig.tight_layout(); fig.savefig(f'{OUT}/fig_desempenho.pdf'); plt.close(fig)

# ---------------------------------------------------- Fig 3: potenciais 1D
fig, ax = plt.subplots(1, 4, figsize=(7.2, 3.0))
xs = np.linspace(-10, 10, 600)
# Morse: De=10, a=0.5, re=3 -> V = 10 (exp(-2a(x-re)) - 2 exp(-a(x-re)))
# DOMINIO REAL: [0,10]  (cf. casos/1D_morse/1D_morse.geo e caso.sif)
xm = np.linspace(0.0, 10.0, 600)
a_morse, re = 0.5, 3.0
V_morse = 10.0*(np.exp(-2*a_morse*(xm-re)) - 2*np.exp(-a_morse*(xm-re)))
ax[0].plot(xm, V_morse, color='0.3', lw=1.4); ax[0].set_ylim(-11.5, 3); ax[0].set_xlim(0, 10)
ax[0].axhline(-8.481033787, color='#c0392b', ls='--', lw=1.2)
ax[0].axhline(-5.818565432, color='#c0392b', ls=':', lw=0.9)
ax[0].set_title('(a) Morse, $a=0.5$, $x_0=3$\n$V_{\\min}=-10$', fontsize=9)
# poco duplo: 0.01 (x^2-4)^2
V_dw = 0.01*(xs**2-4)**2
ax[1].plot(xs, V_dw, color='0.3', lw=1.4); ax[1].set_ylim(-0.05, 1.2); ax[1].set_xlim(-5, 5)
ax[1].axhline(0.224110514, color='#c0392b', ls='--', lw=1.2)
ax[1].axhline(0.600028005, color='#c0392b', ls=':', lw=0.9)
ax[1].set_title('(b) double well\nsplitting $=0.376$', fontsize=9)
# senoidal: A cos(2 pi x / 10), A=0.5
V_sin = 0.5*np.cos(2*np.pi*xs/10.0)
ax[2].plot(xs, V_sin, color='0.3', lw=1.4); ax[2].set_ylim(-0.7, 0.7); ax[2].set_xlim(-10, 10)
ax[2].axhline(-0.213200584, color='#c0392b', ls='--', lw=1.2)
ax[2].set_title('(c) sinusoidal\n$V\\in[-0.5,0.5]$', fontsize=9)
# harmonico
ax[3].plot(xs, xs**2, color='0.3', lw=1.4); ax[3].set_ylim(-1, 12); ax[3].set_xlim(-4, 4)
for E, ls in [(1.000572892, '--'), (3.001197697, ':'), (5.002447062, ':')]:
    ax[3].axhline(E, color='#c0392b', ls=ls, lw=1.2)
ax[3].set_title('(d) harmonic\n$E_n=2n+1$', fontsize=9)
for a in ax:
    a.set_xlabel('$x$', fontsize=9); a.tick_params(labelsize=8)
ax[0].set_ylabel('$V(x)$', fontsize=9)
fig.tight_layout(); fig.savefig(f'{OUT}/fig_potenciais.pdf'); plt.close(fig)

print("  figuras em %s/:" % OUT)
for f in sorted(os.listdir(OUT)):
    print("    %-28s %6.1f KB" % (f, os.path.getsize(os.path.join(OUT, f))/1024))
