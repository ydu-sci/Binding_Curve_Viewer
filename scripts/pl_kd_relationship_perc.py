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
    'axes.labelsize': 'large',
    'axes.titlesize': 'x-large',
    'xtick.labelsize': 'large',
    'ytick.labelsize': 'large',
    'legend.fontsize': 'medium'}
pylab.rcParams.update(params)
color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']


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


# This script generates heatmap of K_d (X), K_i (Y), and equilibrium [PL] (Z)
ki_Y, kd_X= np.mgrid[1:10001, 1:5001]
l0 = 10  # nM
# starting point (pl) 7, 5, 3, 1  nM
# In FP exp., pt should be the total concentration of ligand ([L]0) which equals 10 nM.
# The following equation is from the function `dissociation_constant` in binding_curve_viewer.py
# and is used to calculate the total concentration of ligand to set the configuration of pl.
# lt = (pt + kd - pl) * pl / (pt - pl), substitute the pt with 10 (nM), pl with 7, 5, 3, or 1 (nM)
# lt = (10 + kd - 5) * 5 / (10 - 5)
# lt = kd + 5
p0_7 = (10.0 + kd_X - 7.0) * 7.0 / (10.0 - 7.0)
p0_5 = kd_X + 5.0
p0_3 = (10.0 + kd_X - 3.0) * 3.0 / (10.0 - 3.0)
#p0_1 = (10.0 + kd_X - 1.0) * 1.0 / (10.0 - 1.0)
i0 = 10000.0  # nM, The concentration of the inhibitor is constant (10 uM).
p_eq_7, pl_comp_eq_7, pi_eq_7 = competitive_binding_eq(p0_7, l0, i0, kd_X, ki_Y)
p_eq_5, pl_comp_eq_5, pi_eq_5 = competitive_binding_eq(p0_5, l0, i0, kd_X, ki_Y)
p_eq_3, pl_comp_eq_3, pi_eq_3 = competitive_binding_eq(p0_3, l0, i0, kd_X, ki_Y)
#p_eq_1, pl_comp_eq_1, pi_eq_1 = competitive_binding_eq(p0_1, l0, i0, kd_X, ki_Y)
idx_7 = np.argmin(pl_comp_eq_7, axis=1)
idx_5 = np.argmin(pl_comp_eq_5, axis=1)
idx_3 = np.argmin(pl_comp_eq_3, axis=1)
#idx_1 = np.argmin(pl_comp_eq_1, axis=1)
pl_comp_eq_min_7 = np.min(pl_comp_eq_7, axis=1)
pl_comp_eq_min_5 = np.min(pl_comp_eq_5, axis=1)
pl_comp_eq_min_3 = np.min(pl_comp_eq_3, axis=1)
#pl_comp_eq_min_1 = np.min(pl_comp_eq_1, axis=1)
x_7 = np.take_along_axis(kd_X, np.expand_dims(idx_7, axis=-1), axis=-1).squeeze(axis=-1)
y_7 = np.take_along_axis(ki_Y, np.expand_dims(idx_7, axis=-1), axis=-1).squeeze(axis=-1)
x_5 = np.take_along_axis(kd_X, np.expand_dims(idx_5, axis=-1), axis=-1).squeeze(axis=-1)
y_5 = np.take_along_axis(ki_Y, np.expand_dims(idx_5, axis=-1), axis=-1).squeeze(axis=-1)
x_3 = np.take_along_axis(kd_X, np.expand_dims(idx_3, axis=-1), axis=-1).squeeze(axis=-1)
y_3 = np.take_along_axis(ki_Y, np.expand_dims(idx_3, axis=-1), axis=-1).squeeze(axis=-1)
#x_1 = np.take_along_axis(kd_X, np.expand_dims(idx_1, axis=-1), axis=-1).squeeze(axis=-1)
#y_1 = np.take_along_axis(ki_Y, np.expand_dims(idx_1, axis=-1), axis=-1).squeeze(axis=-1)

fig, ax = plt.subplots(2, 3, figsize=(10, 6))
# Plot contour lines
cs = ax[0, 0].contour(kd_X, ki_Y, pl_comp_eq_3/3.0, colors='k', linewidths=0.5)
cs2 = ax[0, 1].contour(kd_X, ki_Y, pl_comp_eq_5/5.0, colors='k', linewidths=0.5)
cs3 = ax[0, 2].contour(kd_X, ki_Y, pl_comp_eq_7/7.0, colors='k', linewidths=0.5)
# Plot filled contours
cnt = ax[0, 0].contourf(kd_X, ki_Y, pl_comp_eq_3/3.0)
cnt2 = ax[0, 1].contourf(kd_X, ki_Y, pl_comp_eq_5/5.0)
cnt3 = ax[0, 2].contourf(kd_X, ki_Y, pl_comp_eq_7/7.0)
fig.colorbar(cnt, ax=ax[0, 0])
fig.colorbar(cnt2, ax=ax[0, 1])
fig.colorbar(cnt3, ax=ax[0, 2])
ax[0, 0].set_title(r'$[\mathrm{PL}]_0$ = %.1f nM' % 3.0)
ax[0, 0].plot(x_3, y_3, color=color_list[0])
ax[0, 1].set_title(r'$[\mathrm{PL}]_0$ = %.1f nM' % 5.0)
ax[0, 1].plot(x_5, y_5, color=color_list[1])
ax[0, 2].set_title(r'$[\mathrm{PL}]_0$ = %.1f nM' % 7.0)
ax[0, 2].plot(x_7, y_7, color=color_list[2])
ax[1, 0].set_title(r'$[\mathrm{PL}]_0$ = %.1f nM' % 3.0)
ax[1, 0].plot(x_3, pl_comp_eq_min_3/3.0, color=color_list[0])
ax[1, 1].set_title(r'$[\mathrm{PL}]_0$ = %.1f nM' % 5.0)
ax[1, 1].plot(x_5, pl_comp_eq_min_5/5.0, color=color_list[1])
ax[1, 2].set_title(r'$[\mathrm{PL}]_0$ = %.1f nM' % 7.0)
ax[1, 2].plot(x_7, pl_comp_eq_min_7/7.0, color=color_list[2])

ax[0, 0].set_xlabel(r'$K_\mathrm{d}$ ($\mathrm{\mu}$M)')
ax[0, 0].set_ylabel(r'$K_\mathrm{i}$ ($\mathrm{\mu}$M)')
ax[0, 0].set_xticks([0, 2500, 5000])
ax[0, 0].set_xticklabels(['0', '2.5', '5.0'])
ax[0, 0].set_yticks([0, 5000, 10000])
ax[0, 0].set_yticklabels(['0', '5', '10'])
ax[0, 1].set_xlabel(r'$K_\mathrm{d}$ ($\mathrm{\mu}$M)')
ax[0, 1].set_ylabel(r'$K_\mathrm{i}$ ($\mathrm{\mu}$M)')
ax[0, 1].set_xticks([0, 2500, 5000])
ax[0, 1].set_xticklabels(['0', '2.5', '5.0'])
ax[0, 1].set_yticks([0, 5000, 10000])
ax[0, 1].set_yticklabels(['0', '5', '10'])
ax[0, 2].set_xlabel(r'$K_\mathrm{d}$ ($\mathrm{\mu}$M)')
ax[0, 2].set_ylabel(r'$K_\mathrm{i}$ ($\mathrm{\mu}$M)')
ax[0, 2].set_xticks([0, 2500, 5000])
ax[0, 2].set_xticklabels(['0', '2.5', '5.0'])
ax[0, 2].set_yticks([0, 5000, 10000])
ax[0, 2].set_yticklabels(['0', '5', '10'])
ax[1, 0].set_xlabel(r'$K_\mathrm{d}$ (nM)')
ax[1, 0].set_ylabel(r'$[\mathrm{PL}]_\mathrm{eq}/[\mathrm{PL}]_0$')
ax[1, 1].set_xlabel(r'$K_\mathrm{d}$ (nM)')
ax[1, 1].set_ylabel(r'$[\mathrm{PL}]_\mathrm{eq}/[\mathrm{PL}]_0$')
ax[1, 2].set_xlabel(r'$K_\mathrm{d}$ (nM)')
ax[1, 2].set_ylabel(r'$[\mathrm{PL}]_\mathrm{eq}/[\mathrm{PL}]_0$')
ax[1, 2].set_xticks([102, 106, 110])
ax[1, 2].set_xticklabels(['102', '106', '110'])
fig.tight_layout()
fig.savefig('Fig.S11.svg')
plt.close(fig)
