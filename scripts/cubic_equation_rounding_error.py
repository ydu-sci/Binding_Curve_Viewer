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


K1 = 50e-6  # M
c10 = 100  # uM
P0 = 1  # uM

K2 = np.array([0.1, 1, 10, 50]) * 1e-6  # M
C20 = np.arange(0, 101)  # uM
PL1_exact = np.empty([len(K2), len(C20)])

fig, ax = plt.subplots(1, 2, figsize=(10, 5))

for i in range(len(K2)):
    for j in range(len(C20)):
        K2i = K2[i]
        c20 = C20[j]

        a = (K1 + K2i + c10 + c20 - P0)
        b = (K1 * (c20 - P0) + K2i * (c10 - P0) + K1 * K2i)
        g = -K1 * K2i * P0

        t = np.arccos((-2 * a ** 3 + 9 * a * b - 27 * g) / (2 * np.sqrt((a * a - 3 * b) ** 3)))

        P = 2 / 3 * np.sqrt(a * a - 3 * b) * np.cos(t / 3) - a / 3  # [P] / [P]0

        PL1 = c10 * P / (K1 + P)
        # PL2 = c20 * P / (K2i + P)

        PL1_exact[i, j] = PL1

for i in range(len(K2)):
    ax[0].plot(C20, PL1_exact[i, :]/PL1_exact[i, 0], label=r'$K_{d2}$ = %.1f $\mu$M' % (K2[i]*1e6))


# using proper unit
K1 = 50  # uM
c10 = 100  # uM
P0 = 1  # uM

K2 = np.array([0.1, 1, 10, 50])  # uM
C20 = np.arange(0, 101)  # uM
PL1_exact2 = np.empty([len(K2), len(C20)])


for i in range(len(K2)):
    for j in range(len(C20)):
        K2i = K2[i]
        c20 = C20[j]

        a = (K1 + K2i + c10 + c20 - P0)
        b = (K1 * (c20 - P0) + K2i * (c10 - P0) + K1 * K2i)
        g = -K1 * K2i * P0

        t = np.arccos((-2 * a ** 3 + 9 * a * b - 27 * g) / (2 * np.sqrt((a * a - 3 * b) ** 3)))

        P = 2 / 3 * np.sqrt(a * a - 3 * b) * np.cos(t / 3) - a / 3  # [P] / [P]0

        PL1 = c10 * P / (K1 + P)
        # PL2 = c20 * P / (K2i + P)

        PL1_exact2[i, j] = PL1

for i in range(len(K2)):
    ax[1].plot(C20, PL1_exact2[i, :]/PL1_exact2[i, 0], label=r'$K_{d2}$ = %.1f $\mu$M' % K2[i])

ax[0].set_title('Original results')
ax[0].set_xlabel(r'$[L_2]_0$ ($\mu$M)')
ax[0].set_ylabel(r'$[PL_1]_{+L2}/[PL_1]_{-L2}$')
ax[0].legend()
ax[1].set_title('Results calculated by using proper unit')
ax[1].set_xlabel(r'$[L_2]_0$ ($\mu$M)')
ax[1].set_ylabel(r'$[PL_1]_{+L2}/[PL_1]_{-L2}$')
ax[1].legend()
plt.show()