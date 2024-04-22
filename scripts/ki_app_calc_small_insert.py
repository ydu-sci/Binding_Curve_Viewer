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
    'axes.labelsize': 'xx-large',
    'axes.titlesize': 'xx-large',
    'xtick.labelsize': 'xx-large',
    'ytick.labelsize': 'xx-large',}
pylab.rcParams.update(params)


def ic50_curve(p0, l0, pl_kd, pi_kd, resolution=90000):
    a = 1
    b = 0 - (p0 + l0 + pl_kd)
    c = p0 * l0
    pl_eq0 = (0 - b - np.sqrt(b * b - 4 * a * c)) / (2 * a)
    ic50_approx = pi_kd * (1 + l0 / pl_kd)
    plot_ic50_i_log_list = np.linspace(-4, np.log10(ic50_approx * 100), 100)
    plot_ic50_i_list = []
    plot_ic50_i_free_list = []
    plot_ic50_pl_list = []
    for i_log in plot_ic50_i_log_list:
        plot_ic50_i_list.append(10 ** i_log)
    for i in plot_ic50_i_list:
        p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(p0, l0, i, pl_kd, pi_kd)
        plot_ic50_pl_list.append(pl_comp_eq)
        plot_ic50_i_free_list.append(i - pi_eq)
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


def calc_ki_app1(l0, i0, pl_kd):
    # i0 approximates the i_eq, l0 approximates the l_eq
    # the Cheng-Prusoff equation
    # Biochem. Pharmacol. 1973, 22 (23), 3099-3108.
    ki_app1 = i0 / (1 + l0 / pl_kd)
    return ki_app1


pl_kd_list = np.array([50, 5, 50, 5], dtype=float)  # nM
p0_list = np.array([50, 50, 5, 5], dtype=float)  # nM
pi_kd = 100.0  # nM
l0_list = np.arange(100, 501, dtype=float)  # nM


fig, ax = plt.subplots(2, 2, figsize=(10, 10))
for idx in range(len(pl_kd_list)):
    ki1list = []
    for l0 in l0_list:
        plot_ic50_i_list, plot_ic50_i_free_list, plot_ic50_pl_list, apparent_ic50, ic50, pl_eq0 = \
            ic50_curve(p0_list[idx], l0, pl_kd_list[idx], pi_kd)
        pl_eq = pl_eq0 / 2
        ki1 = calc_ki_app1(l0, apparent_ic50, pl_kd_list[idx])
        ki1list.append(ki1)
    ax[0, 0].plot(l0_list, ki1list, label=r'$[P]_0$ = %.1f nM, $K_d$ = %.1f nM' % (p0_list[idx], pl_kd_list[idx]))

fig.tight_layout()
fig.savefig('Figure_5_small_insert.svg')
plt.close(fig)