# MIT License
#
# Copyright (c) 2023 Yu Du
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Binding Curve Viewer: Visualizing the Equilibrium and Kinetics of Protein-Ligand Binding and Competitive Binding"""

__author__ = "Yu Du <orcid.org/0000-0002-4114-396X>"
__version__ = "1.0.0"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

params = {
    'axes.labelsize': 'medium',
    'axes.titlesize': 'medium',
    'xtick.labelsize': 'small',
    'ytick.labelsize': 'small',
    'legend.fontsize': 'small'}
pylab.rcParams.update(params)


def ic50_curve(p0, l0, pl_kd, pi_kd, resolution=90000, left_end=-1):
    a = 1
    b = 0 - (p0 + l0 + pl_kd)
    c = p0 * l0
    pl_eq0 = (0 - b - np.sqrt(b * b - 4 * a * c)) / (2 * a)
    plot_ic50_i_log_list = np.linspace(left_end, 5.5, 100)
    plot_ic50_i_list = []
    plot_ic50_i_free_list = []
    plot_ic50_pl_list = []
    for i_log in plot_ic50_i_log_list:
        plot_ic50_i_list.append(10 ** i_log)
    for i in plot_ic50_i_list:
        p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(p0, l0, i, pl_kd, pi_kd)
        plot_ic50_pl_list.append(pl_comp_eq)
        plot_ic50_i_free_list.append(i - pi_eq)
    ic50_approx = pi_kd * (1 + l0 / pl_kd)
    inner_ic50_i_list = np.linspace(0, 100 * ic50_approx, resolution)
    inner_ic50_i_free_list = []  # this is I free
    inner_ic50_pl_list = []
    for i in inner_ic50_i_list:
        p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(p0, l0, i, pl_kd, pi_kd)
        inner_ic50_pl_list.append(pl_comp_eq)
        inner_ic50_i_free_list.append(i - pi_eq)
    ic50_approx_idx = min(range(len(inner_ic50_pl_list)), key=lambda i: abs(inner_ic50_pl_list[i] - (pl_eq0 / 2)))
    apparent_ic50 = inner_ic50_i_list[ic50_approx_idx]
    ic50 = inner_ic50_i_free_list[ic50_approx_idx]
    return plot_ic50_i_list, plot_ic50_i_free_list, plot_ic50_pl_list, apparent_ic50, ic50, pl_eq0


def competitive_binding_eq(P0, A0, B0, KA, KB):
    # The units of P0, A0, B0, KA and KB are nM
    # The units of  are M
    # Equations from FEBS Lett. 1995, 360 (2), 111-4.
    a = KA + KB + A0 + B0 - P0
    b = KB * (A0 - P0) + KA * (B0 - P0) + KA * KB
    c = -KA * KB * P0
    theta = np.arccos((-2 * a ** 3 + 9 * a * b - 27 * c) / (2 * np.sqrt((a ** 2 - 3 * b) ** 3)))
    p_eq = -a / 3 + (2 / 3) * np.sqrt(a ** 2 - 3 * b) * np.cos(theta / 3)
    pa_eq = p_eq * A0 / (KA + p_eq)
    pb_eq = p_eq * B0 / (KB + p_eq)
    return p_eq, pa_eq, pb_eq


# This script generates the inhibition curves under different condition.
l0 = 10  # nM [L]0 is constant.
# When [P]0=p0_list[i] and Kd=kd_list[i], half ligand (5 nM) is bound.
p0_list = np.array([15, 105, 205, 405], dtype=float)  # nM
kd_list = np.array([10, 100, 200, 400], dtype=float)  # nM
ki_list = np.array([10, 100, 1000, 10000], dtype=float)  # nM

fig, ax = plt.subplots(2, 2, figsize=(10, 10))
color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
for idx in range(len(p0_list)):
    plot_ic50_i_list, plot_ic50_i_free_list, plot_ic50_pl_list, apparent_ic50, ic50, pl_eq0 = \
        ic50_curve(p0_list[idx], l0, kd_list[idx], ki_list[0], left_end=-1)
    ax[0, 0].plot(plot_ic50_i_free_list, plot_ic50_pl_list, color=color_list[idx], linestyle='dotted')
    ax[0, 0].plot(plot_ic50_i_list, plot_ic50_pl_list, color=color_list[idx],
                  label=r'$[\mathrm{P}]_0$ = %.1f, $K_\mathrm{d}$ = %.1f, $\mathrm{IC}_{50}$ = %.1f, App. $\mathrm{IC}_{50}$ = %.1f' %
                        (p0_list[idx], kd_list[idx], ic50, apparent_ic50))
    ax[0, 0].set_xscale('log')

for idx in range(len(p0_list)):
    plot_ic50_i_list, plot_ic50_i_free_list, plot_ic50_pl_list, apparent_ic50, ic50, pl_eq0 = \
        ic50_curve(p0_list[idx], l0, kd_list[idx], ki_list[1], left_end=0)
    ax[0, 1].plot(plot_ic50_i_free_list, plot_ic50_pl_list, color=color_list[idx], linestyle='dotted')
    ax[0, 1].plot(plot_ic50_i_list, plot_ic50_pl_list, color=color_list[idx],
                  label=r'$[\mathrm{P}]_0$ = %.1f, $K_\mathrm{d}$ = %.1f, $\mathrm{IC}_{50}$ = %.1f, App. $\mathrm{IC}_{50}$ = %.1f' %
                        (p0_list[idx], kd_list[idx], ic50, apparent_ic50))
    ax[0, 1].set_xscale('log')

for idx in range(len(p0_list)):
    plot_ic50_i_list, plot_ic50_i_free_list, plot_ic50_pl_list, apparent_ic50, ic50, pl_eq0 = \
        ic50_curve(p0_list[idx], l0, kd_list[idx], ki_list[2], left_end=1)
    ax[1, 0].plot(plot_ic50_i_free_list, plot_ic50_pl_list, color=color_list[idx], linestyle='dotted')
    ax[1, 0].plot(plot_ic50_i_list, plot_ic50_pl_list, color=color_list[idx],
                  label=r'$[\mathrm{P}]_0$ = %.1f, $K_\mathrm{d}$ = %.1f, $\mathrm{IC}_{50}$ = %.1f, App. $\mathrm{IC}_{50}$ = %.1f' %
                        (p0_list[idx], kd_list[idx], ic50, apparent_ic50))
    ax[1, 0].set_xscale('log')

for idx in range(len(p0_list)):
    plot_ic50_i_list, plot_ic50_i_free_list, plot_ic50_pl_list, apparent_ic50, ic50, pl_eq0 = \
        ic50_curve(p0_list[idx], l0, kd_list[idx], ki_list[3], left_end=2)
    ax[1, 1].plot(plot_ic50_i_free_list, plot_ic50_pl_list, color=color_list[idx], linestyle='dotted')
    ax[1, 1].plot(plot_ic50_i_list, plot_ic50_pl_list, color=color_list[idx],
                  label=r'$[\mathrm{P}]_0$ = %.1f, $K_\mathrm{d}$ = %.1f, $\mathrm{IC}_{50}$ = %.1f, App. $\mathrm{IC}_{50}$ = %.1f' %
                        (p0_list[idx], kd_list[idx], ic50, apparent_ic50))
    ax[1, 1].set_xscale('log')

ax[0, 0].set_title(r'$K_\mathrm{i}$ = %.1f nM' % ki_list[0])
ax[0, 0].set_xlabel("Concentration of Free/Total Inhibitor (nM)")
ax[0, 0].set_ylabel("Concentration of Protein-Ligand Complex ([PL], nM)")
ax[0, 0].legend(fontsize='8.2', loc ="upper right")
ax[0, 0].set_ylim([-0.3, 6.8])

ax[0, 1].set_title(r'$K_\mathrm{i}$ = %.1f nM' % ki_list[1])
ax[0, 1].set_xlabel("Concentration of Free/Total Inhibitor (nM)")
ax[0, 1].set_ylabel("Concentration of Protein-Ligand Complex ([PL], nM)")
ax[0, 1].legend(fontsize='8.2', loc ="upper right")
ax[0, 1].set_ylim([-0.3, 6.8])

ax[1, 0].set_title(r'$K_\mathrm{i}$ = %.1f nM' % ki_list[2])
ax[1, 0].set_xlabel("Concentration of Free/Total Inhibitor (nM)")
ax[1, 0].set_ylabel("Concentration of Protein-Ligand Complex ([PL], nM)")
ax[1, 0].legend(fontsize='8.2', loc ="upper right")
ax[1, 0].set_ylim([-0.3, 6.8])

ax[1, 1].set_title(r'$K_\mathrm{i}$ = %.1f nM' % ki_list[3])
ax[1, 1].set_xlabel("Concentration of Free/Total Inhibitor (nM)")
ax[1, 1].set_ylabel("Concentration of Protein-Ligand Complex ([PL], nM)")
ax[1, 1].legend(fontsize='8.2', loc ="upper right")
ax[1, 1].set_ylim([-0.3, 6.8])

fig.tight_layout()
fig.savefig('Fig.S12.png')
plt.close(fig)