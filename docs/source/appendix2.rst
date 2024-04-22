==========
Appendix 2
==========

This appendix includes the background knowledge and equation derivation originally in the supporting information
of the Binding Curve Viewer :ref:`paper<bcv_paper>`.

The calculation of the *K*:sub:`i` and apparent *K*:sub:`i`
===========================================================

In the competitive binding experiment, the binding of the protein and ligand and the binding of the protein and inhibitor are described by the following equations,

.. math::

    \mathrm{P}+\mathrm{L}\overset{K_\mathrm{d}}{\rightleftharpoons}\mathrm{PL}

    \mathrm{P}+\mathrm{I}\overset{K_\mathrm{i}}{\rightleftharpoons}\mathrm{PI}

The :math:`K_\mathrm{d}` is the dissociation constant of the protein and ligand. The :math:`K_\mathrm{i}` is the dissociation
constant of the protein and inhibitor. At equilibrium, the total concentration of the protein (:math:`[\mathrm{P}]_0`)
is the sum of the equilibrium concentration of the protein (:math:`[\mathrm{P}]_\mathrm{eq}`), the equilibrium concentration
of the protein-ligand complex (:math:`[\mathrm{PL}]_\mathrm{eq}`), and the equilibrium concentration of the protein-inhibitor
complex (:math:`[\mathrm{PI}]_\mathrm{eq}`).

.. math:: :label: ap_2.1

    \left[\mathrm{P}\right]_0=\left[\mathrm{P}\right]_\mathrm{eq}+\left[\mathrm{PL}\right]_\mathrm{eq}+\left[\mathrm{PI}\right]_\mathrm{eq}

At equilibrium,

.. math::

    K_\mathrm{d}=\frac{\left[\mathrm{P}\right]_\mathrm{eq}\left[\mathrm{L}\right]_\mathrm{eq}}{\left[\mathrm{PL}\right]_\mathrm{eq}}

.. math:: :label: ap_2.2

    K_\mathrm{i}=\frac{\left[\mathrm{P}\right]_\mathrm{eq}\left[\mathrm{I}\right]_\mathrm{eq}}{\left[\mathrm{PI}\right]_\mathrm{eq}}

Equation 2 can be written as

.. math:: :label: ap_2.3

    \left[\mathrm{PI}\right]_\mathrm{eq}=\frac{\left[\mathrm{P}\right]_\mathrm{eq}\left[\mathrm{I}\right]_\mathrm{eq}}{K_\mathrm{i}}

We substitute equation 3 into equation 1 and obtain

.. math::

    \left[\mathrm{P}\right]_0=\mathrm{P}_\mathrm{eq}+\left[\mathrm{PL}\right]_\mathrm{eq}+\frac{\left[\mathrm{P}\right]_\mathrm{eq}\left[\mathrm{I}\right]_\mathrm{eq}}{K_\mathrm{i}}

.. math:: :label: ap_2.4

    \left[\mathrm{P}\right]_0=\left[\mathrm{P}\right]_\mathrm{eq}\left(1+\frac{\left[\mathrm{I}\right]_\mathrm{eq}}{K_\mathrm{i}}\right)+\left[\mathrm{PL}\right]_\mathrm{eq}

We multiple :math:`\frac{\left[\mathrm{L}\right]_\mathrm{eq}}{\left[\mathrm{PL}\right]_\mathrm{eq}}` on both sides of equation 4 and obtain

.. math::

    \frac{\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_\mathrm{eq}}{\left[\mathrm{PL}\right]_\mathrm{eq}}=\frac{\left[\mathrm{P}\right]_\mathrm{eq}\left[\mathrm{L}\right]_\mathrm{eq}}{\left[\mathrm{PL}\right]_\mathrm{eq}}\left(1+\frac{\left[\mathrm{I}\right]_\mathrm{eq}}{K_\mathrm{i}}\right)+\left[\mathrm{L}\right]_\mathrm{eq}

    \frac{\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_\mathrm{eq}}{\left[\mathrm{PL}\right]_\mathrm{eq}}=K_\mathrm{d}\left(1+\frac{\left[\mathrm{I}\right]_\mathrm{eq}}{K_\mathrm{i}}\right)+\left[\mathrm{L}\right]_\mathrm{eq}

.. math:: :label: ap_2.5

    \left[\mathrm{PL}\right]_\mathrm{eq}=\frac{\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_\mathrm{eq}}{K_\mathrm{d}\left(1+\frac{\left[\mathrm{I}\right]_\mathrm{eq}}{K_\mathrm{i}}\right)+\left[\mathrm{L}\right]_\mathrm{eq}}

Equation 5 can be transformed into equation 6

.. math:: :label: ap_2.6

    K_\mathrm{i}=\frac{\left[\mathrm{I}\right]_\mathrm{eq}}{\frac{\left[\mathrm{L}\right]_\mathrm{eq}\left(\left[\mathrm{P}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}\right)}{\left[\mathrm{PL}\right]_\mathrm{eq}K_\mathrm{d}}-1}

When 50% of the “initial binding” of the protein and ligand is inhibited, the concentration of the free competing inhibitor
(:math:`[\mathrm{I}]_\mathrm{eq}`) is defined as the :math:`\mathrm{IC}_{50}` and the concentration of the total competing inhibitor (:math:`[\mathrm{I}]_0`) is defined as
the apparent :math:`\mathrm{IC}_{50}`. The initial binding means the equilibrium concentration of the protein-ligand complex in
the blank control in the competitive binding of either association or dissociation. In association, the blank
control was the mixture of the protein and the solution with ligand and no inhibitor. In dissociation, the blank
control was the mixture of the equilibrated protein and ligand and the solution without inhibitor. The blank control
is important to eliminate the equilibrium concentration change of the protein-ligand complex upon the volume change
after mixing. At the equilibrium state of 50% inhibition of the protein-ligand binding, the concentration of the
protein-ligand complex (:math:`[\mathrm{PI}]_\mathrm{eq\textnormal{-}50}`) can be written as

.. math:: :label: ap_2.7

    \left[\mathrm{PL}\right]_\mathrm{{eq\textnormal{-}50}}=\frac{\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_\mathrm{{eq\textnormal{-}50}}}{K_\mathrm{d}\left(1+\frac{\left[\mathrm{I}\right]_\mathrm{{eq\textnormal{-}50}}}{K_\mathrm{i}}\right)+\left[\mathrm{L}\right]_\mathrm{{eq\textnormal{-}50}}}

In this work, we used the :ref:`Wang equation <ref1>` to simulate the theoretical inhibition curve, and estimated the
:math:`\mathrm{IC}_{50}` and apparent :math:`\mathrm{IC}_{50}`, without the restrictions of concentrations.

Calculation of the apparent *K*:sub:`i` by the Cheng-Prusoff equation
---------------------------------------------------------------------

At the equilibrium state of 50% inhibition of the protein-ligand binding, if
:math:`\left[\mathrm{L}\right]_0\gg\left[\mathrm{P}\right]_0` and :math:`\left[\mathrm{I}\right]_0\gg\left[\mathrm{P}\right]_0`, or strictly speaking,
:math:`\left[\mathrm{L}\right]_0\gg\left[\mathrm{PL}\right]_\mathrm{{eq\textnormal{-}50}}` and :math:`\left[\mathrm{I}\right]_0\gg\left[\mathrm{PI}\right]_\mathrm{{eq\textnormal{-}50}}`,
:math:`\left[\mathrm{L}\right]_\mathrm{{eq\textnormal{-}50}}=\left[\mathrm{L}\right]_0-\left[\mathrm{PL}\right]_\mathrm{{eq\textnormal{-}50}}\approx\left[\mathrm{L}\right]_0` and
:math:`\mathrm{IC}_{50}=\left[\mathrm{I}\right]_\mathrm{{eq\textnormal{-}50}}=\left[\mathrm{I}\right]_0-\left[\mathrm{PI}\right]_\mathrm{{eq\textnormal{-}50}}\approx\left[\mathrm{I}\right]_0`,
equation 7 approximately equals

.. math::

    \left[\mathrm{PL}\right]_\mathrm{{eq\textnormal{-}50}}\approx\frac{\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_0}{K_\mathrm{d}\left(1+\frac{\left[\mathrm{I}\right]_0}{K_\mathrm{i}}\right)+\left[\mathrm{L}\right]_0}

Thus, :math:`\mathrm{IC}_{50}` is almost equal to :math:`[\mathrm{I}]_0`. From equation 13 in Appendix 1, we know that in the absence of the inhibitor
and when :math:`\left[\mathrm{L}\right]_0\gg\left[\mathrm{P}\right]_0`, the equilibrium concentration of the protein-ligand complex (:math:`[\mathrm{PL}]_\mathrm{{eq\textnormal{-}0}}`) equals

.. math::

    \left[\mathrm{PL}\right]_\mathrm{{eq\textnormal{-}0}}=\frac{\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0}{\left[\mathrm{L}\right]_0+K_\mathrm{d}}

At the time of 50% inhibition of the protein-ligand binding, :math:`\left[\mathrm{PL}\right]_\mathrm{{eq\textnormal{-}0}}=2\times\left[\mathrm{PL}\right]_\mathrm{{eq\textnormal{-}50}}`, so

.. math::

    \frac{\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0}{\left[\mathrm{L}\right]_0+K_\mathrm{d}}=\frac{2\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_0}{K_\mathrm{d}\left(1+\frac{\left[\mathrm{I}\right]_0}{K_\mathrm{i}}\right)+\left[\mathrm{L}\right]_0}

Thus, we obtain the Cheng-Prusoff equation generalized in the competitive binding (equation 8), the :math:`[\mathrm{I}]_0`
in equation 8 is the apparent :math:`\mathrm{IC}_{50}`, which approximately equals the :math:`\mathrm{IC}_{50}` in the condition of
:math:`\left[\mathrm{I}\right]_0\gg\left[\mathrm{PI}\right]_\mathrm{{eq\textnormal{-}50}}`. The apparent :math:`\mathrm{IC}_{50}` can be estimated by the theoretical
inhibition curve.

.. math:: :label: ap_2.8

    K_\mathrm{i}^{\mathrm{app1}}=\frac{\left[\mathrm{I}\right]_0}{1+\frac{\left[\mathrm{L}\right]_0}{K_\mathrm{d}}}

Calculation of the apparent *K*:sub:`i` by the Lin-Riggs equation
-----------------------------------------------------------------

If we substitute equation 5 by :math:`\left[\mathrm{L}\right]_\mathrm{eq}=\left[\mathrm{L}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}`, we obtain

.. math:: :label: ap_2.9

    \left[\mathrm{PL}\right]_\mathrm{eq}=\frac{\left[\mathrm{P}\right]_0\left(\left[\mathrm{L}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}\right)}{K_\mathrm{d}\left(1+\frac{\left[\mathrm{I}\right]_\mathrm{eq}}{K_\mathrm{i}}\right)+\left(\left[\mathrm{L}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}\right)}

We multiple :math:`\frac{1}{\left[\mathrm{L}\right]_0}` on both sides of equation 9 and obtain

.. math:: :label: ap_2.10

    \frac{\left[\mathrm{PL}\right]_\mathrm{eq}}{\left[\mathrm{L}\right]_0}=\frac{\left[\mathrm{P}\right]_0(1-\frac{\left[\mathrm{PL}\right]_\mathrm{eq}}{\left[\mathrm{L}\right]_0})}{K_\mathrm{d}(1+\frac{\left[\mathrm{I}\right]_\mathrm{eq}}{K_\mathrm{i}})+\left[\mathrm{L}\right]_0(1-\frac{\left[\mathrm{PL}\right]_\mathrm{eq}}{\left[\mathrm{L}\right]_0})}

At the equilibrium state of 50% inhibition, equation 10 can be written as

.. math:: :label: ap_2.11

    \frac{\left[\mathrm{PL}\right]_{\mathrm{eq}\textnormal{-}50}}{\left[\mathrm{L}\right]_0}=\frac{\left[\mathrm{P}\right]_0\left(1-\frac{\left[\mathrm{PL}\right]_{\mathrm{eq}\textnormal{-}50}}{\left[\mathrm{L}\right]_0}\right)}{K_\mathrm{d}\left(1+\frac{\left[\mathrm{I}\right]_{\mathrm{eq}\textnormal{-}50}}{K_\mathrm{i}}\right)+\left[\mathrm{L}\right]_0\left(1-\frac{\left[\mathrm{PL}\right]_{\mathrm{eq}\textnormal{-}50}}{\left[\mathrm{L}\right]_0}\right)}

We name the equation (6) in the :ref:`original paper <ref2>` the “Lin-Riggs equation”, which was derived in the
context of the competitive inhibition of the protein and labelled nucleic acid. In contrast to
:math:`\frac{\left[\mathrm{PL}\right]_\mathrm{{eq\textnormal{-}50}}}{\left[\mathrm{PL}\right]_\mathrm{{eq\textnormal{-}0}}}=0.5` in the previous section,
:math:`\frac{\left[\mathrm{PL}\right]_\mathrm{{eq\textnormal{-}50}}}{\left[\mathrm{L}\right]_0}=0.5` in the competitive inhibition of the
protein and labelled nucleic acid. Thus, equation 11 can be converted to

.. math:: :label: ap_2.12

    K_\mathrm{i}=\frac{2\left[\mathrm{I}\right]_\mathrm{{eq\textnormal{-}50}}K_\mathrm{d}}{2\left[\mathrm{P}\right]_0-\left[\mathrm{L}\right]_0-2K_\mathrm{d}}

The usage condition of the Lin-Riggs equation is :math:`\left[\mathrm{L}\right]_0\approx\left[\mathrm{PL}\right]_\mathrm{eq}`
and :math:`\left[\mathrm{I}\right]_0\gg\left[\mathrm{PI}\right]_\mathrm{eq}`. In this condition,
:math:`\mathrm{IC}_{50}=\left[\mathrm{I}\right]_\mathrm{{eq\textnormal{-}50}}=\left[\mathrm{I}\right]_0-\left[\mathrm{PI}\right]_\mathrm{{eq\textnormal{-}50}}\approx\left[\mathrm{I}\right]_0`,
we can thus use equation 13 to calculate the apparent :math:`K_\mathrm{i}`,

.. math:: :label: ap_2.13

    K_\mathrm{i}^{\mathrm{app2}}=\frac{2\left[\mathrm{I}\right]_0K_\mathrm{d}}{2\left[\mathrm{P}\right]_0-\left[\mathrm{L}\right]_0-2K_\mathrm{d}}

The :math:`[\mathrm{I}]_0` in equation 13 is the apparent :math:`\mathrm{IC}_{50}`, which approximately equals the :math:`\mathrm{IC}_{50}`. Thus,
equation 13 is essentially the same as the Lin-Riggs equation. The apparent :math:`\mathrm{IC}_{50}` can be
estimated by the theoretical inhibition curve.

It is worth noting that the experimental conditions were different between the calculations
using the Cheng-Prusoff equation and the Lin-Riggs equation, because the :math:`[\mathrm{PL}]_\mathrm{{eq\textnormal{-}0}}` should
be greater than the half of the :math:`[\mathrm{L}]_0` in the calculations using the Lin-Riggs equation.
Namely, in the conditions of :math:`[\mathrm{P}]_0 = 50\ \mathrm{nM}` and :math:`K_\mathrm{d} = 50\ \mathrm{nM}`,
:math:`[\mathrm{P}]_0 = 5\ \mathrm{nM}` and :math:`K_\mathrm{d} = 50\ \mathrm{nM}`,
and :math:`[\mathrm{P}]_0 = 5\ \mathrm{nM}` and :math:`K_\mathrm{d} = 5\ \mathrm{nM}`,
the :math:`[\mathrm{PL}]_\mathrm{{eq\textnormal{-}0}}` was smaller than :math:`\frac{\left[\mathrm{L}\right]_0}{2}`.
When generating the theoretical inhibition curve with increasing :math:`[\mathrm{L}]_0`, the range of
the :math:`[\mathrm{L}]_0` was also restricted, because the :math:`[\mathrm{PL}]_\mathrm{{eq\textnormal{-}0}}`
should be greater than the half of the :math:`[\mathrm{L}]_0`.

Calculation of the *K*:sub:`i` by the Wang’s group equation
-----------------------------------------------------------

We used the :ref:`Wang’s group equation <ref3>` (equation 14) to calculate the :math:`K_\mathrm{i}`.
At the equilibrium state of 50% inhibition, :math:`\frac{\left[\mathrm{PL}\right]_\mathrm{{eq\textnormal{-}50}}}{\left[\mathrm{PL}\right]_\mathrm{{eq\textnormal{-}0}}}=0.5`.
The apparent :math:`\mathrm{IC}_{50}` can be estimated by the theoretical inhibition curve.

.. math:: :label: ap_2.14

    K_\mathrm{i}=\frac{\left[\mathrm{I}\right]_\mathrm{{eq\textnormal{-}50}}}{\frac{\left[\mathrm{L}\right]_\mathrm{{eq\textnormal{-}50}}}{K_\mathrm{d}}+\frac{\left[\mathrm{P}\right]_\mathrm{{eq\textnormal{-}0}}}{K_\mathrm{d}}+1}

For the details of the calculation, see the :ref:`original paper <ref3>` or the function calc_ki in
our script ki_app_calc.py

