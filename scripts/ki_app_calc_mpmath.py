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
from mpmath import mp

mp.dps = 50


def ic50_curve(p0, l0, pl_kd, pi_kd, resolution=90000):
    a = mp.mpf(1)
    b = mp.fneg(mp.fsum([p0, l0, pl_kd]))
    c = mp.fmul(p0, l0)
    pl_eq0 = mp.fdiv(mp.fneg(mp.fadd(b, mp.sqrt(mp.fsub(mp.fmul(b, b),  mp.fmul(mp.fmul(4, a), c))))), mp.fmul(2, a))
    ic50_approx = mp.fmul(pi_kd, mp.fadd(1, mp.fdiv(l0, pl_kd)))
    inner_ic50_i_list = mp.linspace(0, mp.fmul(100, ic50_approx), resolution)
    inner_ic50_i_free_list = []  # this is I free
    inner_ic50_pl_list = []
    for i in inner_ic50_i_list:
        p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(p0, l0, i, pl_kd, pi_kd)
        inner_ic50_pl_list.append(pl_comp_eq)
        inner_ic50_i_free_list.append(mp.fsub(i, pi_eq))
    ic50_approx_idx = min(range(len(inner_ic50_pl_list)), key=lambda i: mp.fabs(mp.fsub(inner_ic50_pl_list[i],
                                                                                        mp.fdiv(pl_eq0, 2))))
    apparent_ic50 = inner_ic50_i_list[ic50_approx_idx]
    ic50 = inner_ic50_i_free_list[ic50_approx_idx]
    return apparent_ic50, ic50, pl_eq0


def competitive_binding_eq(P0, A0, B0, KA, KB):
    # The units of P0, A0, B0, KA and KB are nM
    # The units of  are M
    # Equations from FEBS Lett. 1995, 360 (2), 111-4.
    a = mp.fsub(mp.fsum([KA, KB, A0, B0]), P0)
    b = mp.fsum([mp.fmul(KB, mp.fsub(A0, P0)), mp.fmul(KA, mp.fsub(B0, P0)), mp.fmul(KA, KB)])
    c = mp.fneg(mp.fmul(KA, mp.fmul(KB, P0)))
    theta = mp.acos(mp.fdiv(mp.fadd(mp.fneg(mp.fmul(2, mp.power(a, 3))), mp.fsub(mp.fmul(9, mp.fmul(a, b)), mp.fmul(27, c))),
            mp.fmul(2, mp.sqrt(mp.power(mp.fsub(mp.fmul(a, a), mp.fmul(3, b)), 3)))))
    p_eq = mp.fadd(mp.fdiv(mp.fneg(a), 3),
                   mp.fmul(mp.fdiv(2, 3), mp.fmul(mp.sqrt(mp.fsub(mp.fmul(a, a), mp.fmul(3, b))), mp.cos(mp.fdiv(theta, 3)))))
    pa_eq = mp.fdiv(mp.fmul(p_eq, A0), mp.fadd(KA, p_eq))
    pb_eq = mp.fdiv(mp.fmul(p_eq, B0), mp.fadd(KB, p_eq))
    return p_eq, pa_eq, pb_eq


def calc_ki(pt, lt, it, pl_kd, pl_eq0):
    # use the solutions in Anal. Biochem. 2004, 332 (2), 261-73.
    # pt - pl_1 is the P0 in original paper.
    # pl_1 is the PL0 in original paper.
    P0 = mp.fsub(pt, pl_eq0)
    PL50 = mp.fdiv(pl_eq0, 2)
    L50 = mp.fsub(lt, PL50)
    I50 = mp.fsum([mp.fsub(it, pt), mp.fmul(pl_kd, mp.fdiv(PL50, L50)), PL50])
    ki = mp.fdiv(I50, mp.fsum([mp.fdiv(L50, pl_kd), mp.fdiv(P0, pl_kd), 1]))
    return ki


#pl_kd_list = [mp.mpf(50), mp.mpf(5), mp.mpf(50), mp.mpf(5)]  # nM
pl_kd_list = [mp.mpf(50)]  # nM
#p0_list = [mp.mpf(50), mp.mpf(50), mp.mpf(5), mp.mpf(50)]  # nM
p0_list = [mp.mpf(50)]  # nM
#p0_list_plot = np.array([50, 50, 5, 5], dtype=float)  # nM
p0_list_plot = np.array([50], dtype=float)  # nM
#pl_kd_list_plot = np.array([50, 5, 50, 5], dtype=float)  # nM
pl_kd_list_plot = np.array([50], dtype=float)  # nM
pi_kd = mp.mpf(100)  # nM
l0_list = mp.arange(1, 501)
l0_list_plot = np.arange(1, 501)


fig, ax = plt.subplots(1, 2, figsize=(10, 5))
for idx in range(len(pl_kd_list)):
    ki3list = []
    for l0 in l0_list:
        apparent_ic50, ic50, pl_eq0 = ic50_curve(p0_list[idx], l0, pl_kd_list[idx], pi_kd)
        ki3 = calc_ki(p0_list[idx], l0, apparent_ic50, pl_kd_list[idx], pl_eq0)
        ki3list.append(ki3)
    ki3list_plot = np.array([mp.nstr(i, 6) for i in ki3list], dtype=float)
    ax[0].plot(l0_list_plot, ki3list_plot, label=r'$[{\mathrm{P}}]_0$ = %.1f nM, $K_{\mathrm{d}}$ = %.1f nM' % (p0_list_plot[idx], pl_kd_list_plot[idx]))

for idx in range(len(pl_kd_list)):
    ki3list = []
    for l0 in l0_list:
        apparent_ic50, ic50, pl_eq0 = ic50_curve(p0_list[idx], l0, pl_kd_list[idx], pi_kd, 180000)
        ki3 = calc_ki(p0_list[idx], l0, apparent_ic50, pl_kd_list[idx], pl_eq0)
        ki3list.append(ki3)
    ki3list_plot = np.array([mp.nstr(i, 6) for i in ki3list], dtype=float)
    ax[1].plot(l0_list_plot, ki3list_plot, label=r'$[{\mathrm{P}}]_0$ = %.1f nM, $K_{\mathrm{d}}$ = %.1f nM' % (p0_list_plot[idx], pl_kd_list_plot[idx]))

ax[0].set_title("Wang's group equation")
ax[0].set_xlabel(r'$[{\mathrm{L}}]_{\mathrm{0}}$ (nM)')
ax[0].set_ylabel(r'$K_{\mathrm{i}}$ (nM)')
ax[0].legend()

ax[1].set_title(r"Wang's group equation with more accurate apparent ${\mathrm{IC}}_{50}$")
ax[1].set_xlabel(r'$[{\mathrm{L}}]_{\mathrm{0}}$ (nM)')
ax[1].set_ylabel(r'$K_{\mathrm{i}}$ (nM)')
ax[1].legend()

fig.tight_layout()
fig.savefig('Fig.S10-mpmath.svg')
plt.close(fig)