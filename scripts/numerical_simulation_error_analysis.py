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


def calc_disso_t99(koff):
    disso_t99 = np.log(100) / koff
    return disso_t99


def second_order_kinetics_t99(p0, l0, koff, kon):
    kon = kon * 1e-9  # convert the unit from M^-1 s^-1 to nM^-1 s^-1
    a = 1
    b = 0 - (p0 + l0 + koff / kon)
    c = p0 * l0
    pl_1 = (0 - b - np.sqrt(b * b - 4 * a * c)) / (2 * a)  # pl_1 is [PL]eq, pl_1 < pl_2
    pl_2 = (0 - b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
    A = 1 / (pl_1 - pl_2)
    num_array = np.arange(0, 100, 1) / 100
    asso_list = pl_1 * num_array
    time_list = A / kon * np.log((pl_2 * (pl_1 - asso_list)) / (pl_1 * (pl_2 - asso_list)))
    t99 = time_list[99]
    return t99, pl_1


def competitive_binding_eq(P0, A0, B0, KA, KB):
    # The units of P0, A0, and B0 are nM
    # The units of KA and KB are M
    # Equations from FEBS Lett. 1995, 360 (2), 111-4.
    # KA = koff1 / kon1
    # KB = koff2 / kon2
    KA = KA * 1e9  # Change unit of M to nM
    KB = KB * 1e9
    a = KA + KB + A0 + B0 - P0
    b = KB * (A0 - P0) + KA * (B0 - P0) + KA * KB
    c = -KA * KB * P0
    theta = np.arccos((-2 * a ** 3 + 9 * a * b - 27 * c) / (2 * np.sqrt((a ** 2 - 3 * b) ** 3)))
    p_eq = -a / 3 + (2 / 3) * np.sqrt(a ** 2 - 3 * b) * np.cos(theta / 3)
    pa_eq = p_eq * A0 / (KA + p_eq)
    pb_eq = p_eq * B0 / (KB + p_eq)
    return p_eq, pa_eq, pb_eq


def diff_calc_association_eq(p0, l0, i0, pl_koff, pl_kon, pi_koff, pi_kon,
                             t99_short, t99_long, resolution=200000, t99Ntimes=10):
    # The units of p0, l0, and i0 are nM
    # The units of pl_koff and pi_koff are s-1
    # The units of pl_kon and pi_kon are M-1 s-1
    pl_kon = pl_kon * 1e-9
    pi_kon = pi_kon * 1e-9
    t99_2 = t99_long * t99Ntimes
    time_loop1, dt1 = np.linspace(0, t99_short, resolution, retstep=True)
    time_loop2, dt2 = np.linspace(t99_short, t99_2, resolution, retstep=True)
    con_pl_0 = p0 * l0 * pl_kon * dt1
    coff_pl_0 = con_pl_0 * pl_koff * dt1
    con_pi_0 = p0 * i0 * pi_kon * dt1
    coff_pi_0 = con_pi_0 * pi_koff * dt1
    con0 = con_pl_0 + con_pi_0
    coff0 = coff_pl_0 + coff_pi_0
    conlist_p, cofflist_p, consumlist_p, coffsumlist_p = [con0], [coff0], [con0], [coff0]
    conlist_l, cofflist_l, consumlist_l, coffsumlist_l = [con_pl_0], [coff_pl_0], [con_pl_0], [coff_pl_0]
    conlist_i, cofflist_i, consumlist_i, coffsumlist_i = [con_pi_0], [coff_pi_0], [con_pi_0], [coff_pi_0]
    for s in range(2, resolution):
        con_l = (p0 - consumlist_p[-1] + coffsumlist_p[-1]) * \
                (l0 - consumlist_l[-1] + coffsumlist_l[-1]) * pl_kon * dt1
        con_i = (p0 - consumlist_p[-1] + coffsumlist_p[-1]) * \
                (i0 - consumlist_i[-1] + coffsumlist_i[-1]) * pi_kon * dt1
        con_p = con_l + con_i
        conlist_p.append(con_p)
        conlist_l.append(con_l)
        conlist_i.append(con_i)
        consumlist_p.append(consumlist_p[-1] + con_p)
        consumlist_l.append(consumlist_l[-1] + con_l)
        consumlist_i.append(consumlist_i[-1] + con_i)
        coff_l = (consumlist_l[-1] - coffsumlist_l[-1]) * pl_koff * dt1
        coff_i = (consumlist_i[-1] - coffsumlist_i[-1]) * pi_koff * dt1
        coff_p = coff_l + coff_i
        cofflist_p.append(coff_p)
        cofflist_l.append(coff_l)
        cofflist_i.append(coff_i)
        coffsumlist_p.append(coffsumlist_p[-1] + coff_p)
        coffsumlist_l.append(coffsumlist_l[-1] + coff_l)
        coffsumlist_i.append(coffsumlist_i[-1] + coff_i)
    for s in range(1, resolution):
        con_l = (p0 - consumlist_p[-1] + coffsumlist_p[-1]) * \
                (l0 - consumlist_l[-1] + coffsumlist_l[-1]) * pl_kon * dt2
        con_i = (p0 - consumlist_p[-1] + coffsumlist_p[-1]) * \
                (i0 - consumlist_i[-1] + coffsumlist_i[-1]) * pi_kon * dt2
        con_p = con_l + con_i
        conlist_p.append(con_p)
        conlist_l.append(con_l)
        conlist_i.append(con_i)
        consumlist_p.append(consumlist_p[-1] + con_p)
        consumlist_l.append(consumlist_l[-1] + con_l)
        consumlist_i.append(consumlist_i[-1] + con_i)
        coff_l = (consumlist_l[-1] - coffsumlist_l[-1]) * pl_koff * dt2
        coff_i = (consumlist_i[-1] - coffsumlist_i[-1]) * pi_koff * dt2
        coff_p = coff_l + coff_i
        cofflist_p.append(coff_p)
        cofflist_l.append(coff_l)
        cofflist_i.append(coff_i)
        coffsumlist_p.append(coffsumlist_p[-1] + coff_p)
        coffsumlist_l.append(coffsumlist_l[-1] + coff_l)
        coffsumlist_i.append(coffsumlist_i[-1] + coff_i)
    pl_list = [0]  # length should be 2 * resolution - 1
    # the length of consumlist_l is 2 * resolution - 2
    for i in range(0, 2 * resolution - 2):
        pl_list.append(consumlist_l[i] - coffsumlist_l[i])
    return pl_list[-1]

def diff_calc_dissociation_eq(pl0, p0, l0, i0, pl_koff, pl_kon, pi_koff, pi_kon,
                              t99_short, t99_long, resolution=400000, t99Ntimes=10):
    # The units of p0, l0, and i0 are nM
    # The units of pl_koff and pi_koff are s-1
    # The units of pl_kon and pi_kon are M-1 s-1
    pi0 = 0
    pl_kon = pl_kon * 1e-9
    pi_kon = pi_kon * 1e-9
    t99_2 = t99_long * t99Ntimes
    time_loop1, dt1 = np.linspace(0, t99_short, resolution, retstep=True)
    time_loop2, dt2 = np.linspace(t99_short, t99_2, resolution, retstep=True)
    con_pl_0 = p0 * l0 * pl_kon * dt1
    coff_pl_0 = (con_pl_0 + pl0) * pl_koff * dt1
    con_pi_0 = p0 * i0 * pi_kon * dt1
    coff_pi_0 = (con_pi_0 + pi0) * pi_koff * dt1
    con0 = con_pl_0 + con_pi_0
    coff0 = coff_pl_0 + coff_pi_0
    conlist_p, cofflist_p, consumlist_p, coffsumlist_p = [con0], [coff0], [con0], [coff0]
    conlist_l, cofflist_l, consumlist_l, coffsumlist_l = [con_pl_0], [coff_pl_0], [con_pl_0], [coff_pl_0]
    conlist_i, cofflist_i, consumlist_i, coffsumlist_i = [con_pi_0], [coff_pi_0], [con_pi_0], [coff_pi_0]
    for s in range(2, resolution):
        con_l = (p0 - consumlist_p[-1] + coffsumlist_p[-1]) * \
                (l0 - consumlist_l[-1] + coffsumlist_l[-1]) * pl_kon * dt1
        con_i = (p0 - consumlist_p[-1] + coffsumlist_p[-1]) * \
                (i0 - consumlist_i[-1] + coffsumlist_i[-1]) * pi_kon * dt1
        con_p = con_l + con_i
        conlist_p.append(con_p)
        conlist_l.append(con_l)
        conlist_i.append(con_i)
        consumlist_p.append(consumlist_p[-1] + con_p)
        consumlist_l.append(consumlist_l[-1] + con_l)
        consumlist_i.append(consumlist_i[-1] + con_i)
        coff_l = (pl0 + consumlist_l[-1] - coffsumlist_l[-1]) * pl_koff * dt1
        coff_i = (pi0 + consumlist_i[-1] - coffsumlist_i[-1]) * pi_koff * dt1
        coff_p = coff_l + coff_i
        cofflist_p.append(coff_p)
        cofflist_l.append(coff_l)
        cofflist_i.append(coff_i)
        coffsumlist_p.append(coffsumlist_p[-1] + coff_p)
        coffsumlist_l.append(coffsumlist_l[-1] + coff_l)
        coffsumlist_i.append(coffsumlist_i[-1] + coff_i)
    for s in range(1, resolution):
        con_l = (p0 - consumlist_p[-1] + coffsumlist_p[-1]) * \
                (l0 - consumlist_l[-1] + coffsumlist_l[-1]) * pl_kon * dt2
        con_i = (p0 - consumlist_p[-1] + coffsumlist_p[-1]) * \
                (i0 - consumlist_i[-1] + coffsumlist_i[-1]) * pi_kon * dt2
        con_p = con_l + con_i
        conlist_p.append(con_p)
        conlist_l.append(con_l)
        conlist_i.append(con_i)
        consumlist_p.append(consumlist_p[-1] + con_p)
        consumlist_l.append(consumlist_l[-1] + con_l)
        consumlist_i.append(consumlist_i[-1] + con_i)
        coff_l = (pl0 + consumlist_l[-1] - coffsumlist_l[-1]) * pl_koff * dt2
        coff_i = (pi0 + consumlist_i[-1] - coffsumlist_i[-1]) * pi_koff * dt2
        coff_p = coff_l + coff_i
        cofflist_p.append(coff_p)
        cofflist_l.append(coff_l)
        cofflist_i.append(coff_i)
        coffsumlist_p.append(coffsumlist_p[-1] + coff_p)
        coffsumlist_l.append(coffsumlist_l[-1] + coff_l)
        coffsumlist_i.append(coffsumlist_i[-1] + coff_i)
    pl_list = [pl0]  # length should be 2 * resolution - 1
    # the length of consumlist_l is 2 * resolution - 2
    for i in range(0, 2 * resolution - 2):
        pl_list.append(pl0 + consumlist_l[i] - coffsumlist_l[i])
    return pl_list[-1]

p0 = [0.1, 1, 10, 100]  # nM
l0 = [1, 10, 100, 1000]  # nM
i0 = [10, 100, 1000, 10000]  # nM
pl_kon = [1e7, 1e6, 1e7, 1e6]  # M-1 s-1
pl_koff = [1e-2, 1e-2, 1e-1, 1e-1]  # s-1
pi_kon = [1e6, 1e7, 1e6, 1e6]  # M-1 s-1
pi_koff = [1e-2, 1e-1, 1e-1, 1]  # s-1
# association analysis
for i in range(4):
    pl_kd = pl_koff[i] / pl_kon[i]
    pi_kd = pi_koff[i] / pi_kon[i]
    t99_pl, pl_eq_second_order = second_order_kinetics_t99(p0[i], l0[i], pl_koff[i], pl_kon[i])
    t99_pi, pi_eq_second_order = second_order_kinetics_t99(p0[i], i0[i], pi_koff[i], pi_kon[i])
    t99_short = None
    t99_long = None
    if t99_pl >= t99_pi:
        t99_long = t99_pl
        t99_short = t99_pi
    else:
        t99_long = t99_pi
        t99_short = t99_pl

    p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(p0[i], l0[i], i0[i], pl_kd, pi_kd)
    diff_eq = diff_calc_association_eq(p0[i], l0[i], i0[i], pl_koff[i], pl_kon[i], pi_koff[i], pi_kon[i],
                                       t99_short, t99_long)
    if abs(diff_eq - pl_comp_eq) > 0.001 * pl_comp_eq:
        print('p0, l0, i0, pl_koff, pl_kon, pi_koff, pi_kon')
        print(p0[i], l0[i], i0[i], pl_koff[i], pl_kon[i], pi_koff[i], pi_kon[i])
        print(diff_eq, pl_comp_eq)

p0_disso = [0.2, 2, 20, 200]  # nM
l0_disso = [2, 20, 200, 2000]  # nM
i0_disso = [20, 200, 2000, 20000]  # nM
volume_ratio = 1
# dissociation analysis
for i in range(4):
    pl_kd = pl_koff[i] / pl_kon[i]
    pi_kd = pi_koff[i] / pi_kon[i]
    t99_pl, pl_eq_second_order = second_order_kinetics_t99(p0_disso[i], l0_disso[i], pl_koff[i], pl_kon[i])
    t99_pi, pi_eq_second_order = second_order_kinetics_t99(p0_disso[i], i0_disso[i], pi_koff[i], pi_kon[i])
    pl_disso_t99 = calc_disso_t99(pl_koff[i])
    t99_short = None
    t99_long = None
    if pl_disso_t99 >= t99_pi:
        t99_long = pl_disso_t99
        t99_short = t99_pi
    else:
        t99_long = t99_pi
        t99_short = pl_disso_t99
    p_total_mix = p0_disso[i] * volume_ratio / (1 + volume_ratio)
    l_total_mix = l0_disso[i] * volume_ratio / (1 + volume_ratio)
    i_total_mix = i0_disso[i] / (1 + volume_ratio)
    p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(p_total_mix, l_total_mix, i_total_mix, pl_kd, pi_kd)

    pl_mix = pl_eq_second_order * volume_ratio / (1 + volume_ratio)
    p_mix = (p0_disso[i] - pl_eq_second_order) * volume_ratio / (1 + volume_ratio)
    l_mix = (l0_disso[i] - pl_eq_second_order) * volume_ratio / (1 + volume_ratio)
    i_mix = i0_disso[i] / (1 + volume_ratio)
    diff_eq_disso = diff_calc_dissociation_eq(pl_mix, p_mix, l_mix, i_mix, pl_koff[i], pl_kon[i], pi_koff[i],
                                              pi_kon[i], t99_short, t99_long)
    if abs(diff_eq_disso - pl_comp_eq) > 0.001 * pl_comp_eq:
        print('p0_disso, l0_disso, i0_disso, pl_koff, pl_kon, pi_koff, pi_kon, volume_ratio')
        print(p0_disso[i], l0_disso[i], i0_disso[i], pl_koff[i], pl_kon[i], pi_koff[i], pi_kon[i], volume_ratio)
        print(diff_eq_disso, pl_comp_eq)

