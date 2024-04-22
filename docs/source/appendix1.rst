==========
Appendix 1
==========

This appendix includes the background knowledge and equation derivation originally in the supporting information
of the Binding Curve Viewer :ref:`paper<bcv_paper>`.

Solutions of differential equations and equilibrium concentrations
==================================================================

Throughout this section, we describe the binding process initiated by mixing an arbitrary volume of protein with
an arbitrary volume of ligand. This meant that there was no protein-ligand complex at the start of the binding
process. The binding of the protein and ligand was 1:1 stoichiometry and had no cooperativity. In the dissociation
process, we assumed that there was no rebinding.

.. _kinetics_2nd:

Kinetics of the second-order binding process
--------------------------------------------

The binding of the protein and ligand is described by the following equation,

.. math::

    \mathrm{P} + \mathrm{L} \mathrel{\mathop{\rightleftarrows}^{k_{\mathrm{off}}}_{k_{\mathrm{on}}}} \mathrm{PL}

The concentration of the protein-ligand complex changes during the binding process and

.. math:: :label: ap_1.1

    \frac{d\left[\mathrm{PL}\right]}{dt}=k_{\mathrm{on}}\left[\mathrm{P}\right]\left[\mathrm{L}\right]-k_{\mathrm{off}}\left[\mathrm{PL}\right]

where the :math:`[\mathrm{P}]` and :math:`[\mathrm{L}]` are the concentrations of the unbound protein and ligand,
the :math:`k_{\mathrm{on}}` and :math:`k_{\mathrm{off}}` is respectively the association rate constant and dissociation rate constant. The concentration
of the total protein and ligand, :math:`[\mathrm{P}]_0` and :math:`[\mathrm{L}]_0`, and the concentration of the
unbound protein and ligand, :math:`[\mathrm{P}]` and :math:`[\mathrm{L}]` are related by

.. math::

    \left[\mathrm{P}\right]=\left[\mathrm{P}\right]_0-\left[\mathrm{PL}\right]

    \left[\mathrm{L}\right]=\left[\mathrm{L}\right]_0-\left[\mathrm{PL}\right]

Thus,

.. math::

    \frac{d\left[\mathrm{PL}\right]}{dt}=k_{\mathrm{on}}\left(\left[\mathrm{P}\right]_0-\left[\mathrm{\mathrm{PL}}\right]\right)\left(\left[\mathrm{L}\right]_0-\left[\mathrm{PL}\right]\right)-k_{\mathrm{off}}\left[\mathrm{PL}\right]

    \frac{d\left[\mathrm{PL}\right]}{dt}=k_{\mathrm{on}}\left[\mathrm{PL}\right]^2-\left(k_{\mathrm{on}}\left[\mathrm{P}\right]_0+k_{\mathrm{on}}\left[\mathrm{L}\right]_0+k_{\mathrm{off}}\right)\left[\mathrm{PL}\right]+k_{\mathrm{on}}\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_0


.. math:: :label: ap_1.2

    \frac{d\left[\mathrm{PL}\right]}{k_{\mathrm{on}}dt}=\left[\mathrm{PL}\right]^2-\left(\left[\mathrm{P}\right]_0+\left[L\right]_0+\frac{k_{\mathrm{off}}}{k_{\mathrm{on}}}\right)\left[\mathrm{PL}\right]+\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_0

The coefficients and discriminant of the right-hand side quadratic equation are

.. math::

    a=1,\ b=-\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0+\frac{k_{\mathrm{off}}}{k_{\mathrm{on}}}\right)=-\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0+K_\mathrm{d}\right),\ c=\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_0

    \Delta=b^2-4ac=\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0+K_\mathrm{d}\right)^2-4\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_0=\left(\left[\mathrm{P}\right]_0-\left[\mathrm{L}\right]_0\right)^2+{K_\mathrm{d}}^2+2\left[\mathrm{P}\right]_0K_\mathrm{d}+2\left[\mathrm{L}\right]_0K_\mathrm{d}>0

Thus, two distinct real roots are

.. math::

    \left[\mathrm{PL}\right]_1=\frac{-b-\sqrt{b^2-4ac}}{2a}=\frac{\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0+K_\mathrm{d}\right)-\sqrt{\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0+K_\mathrm{d}\right)^2-4\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_0}}{2}

    \left[\mathrm{PL}\right]_2=\frac{-b+\sqrt{b^2-4ac}}{2a}

Then equation 2 may be written as

.. math::

    \frac{d\left[\mathrm{PL}\right]}{k_{\mathrm{on}}dt}=(\left[\mathrm{PL}\right]-\left[\mathrm{PL}\right]_1)(\left[\mathrm{PL}\right]-\left[\mathrm{PL}\right]_2)

    \frac{d\left[\mathrm{PL}\right]}{\left(\left[\mathrm{PL}\right]-\left[\mathrm{PL}\right]_1\right)\left(\left[\mathrm{PL}\right]-\left[\mathrm{PL}\right]_2\right)}=k_{\mathrm{on}}dt

    \frac{1}{\left[\mathrm{PL}\right]_1-\left[\mathrm{PL}\right]_2}(\frac{1}{\left[\mathrm{PL}\right]-\left[\mathrm{PL}\right]_1}-\frac{1}{\left[\mathrm{PL}\right]-\left[\mathrm{PL}\right]_2})d[\mathrm{PL}]=k_{\mathrm{on}}dt

    \frac{d\left[\mathrm{PL}\right]}{\left[\mathrm{PL}\right]-\left[\mathrm{PL}\right]_1}-\frac{d\left[\mathrm{PL}\right]}{\left[\mathrm{PL}\right]-\left[\mathrm{PL}\right]_2}=k_{\mathrm{on}}\left(\left[\mathrm{PL}\right]_1-\left[\mathrm{PL}\right]_2\right)dt

    \ln{\left|\left[\mathrm{PL}\right]-\left[\mathrm{PL}\right]_1\right|}+\mathrm{C}_1-\ln{\left|\left[\mathrm{PL}\right]-\left[\mathrm{PL}\right]_2\right|}+\mathrm{C}_2=k_{\mathrm{on}}\left(\left[\mathrm{PL}\right]_1-\left[\mathrm{PL}\right]_2\right)t+\mathrm{C}_3

Because :math:`b<0`, :math:`[\mathrm{PL}]_1<[\mathrm{PL}]_2`. From the expression of :math:`[\mathrm{PL}]_1`, we know that :math:`[\mathrm{PL}]_1>0`.
When :math:`[\mathrm{PL}]` increases from 0 to :math:`[\mathrm{PL}]_1`,  :math:`\frac{d\left[\mathrm{PL}\right]}{dt}` decreases from the maximum value to 0, i.e., :math:`\frac{d\left[\mathrm{PL}\right]}{dt}\rightarrow0`, the binding process approaches equilibrium, so
:math:`\left[\mathrm{PL}\right]\in\left[0,\ \left[\mathrm{PL}\right]_1\right)`. Thus, the above equation gives

.. math::

    \ln{\left(\left[\mathrm{PL}\right]_1-\left[\mathrm{PL}\right]\right)}-\ln{\left(\left[\mathrm{PL}\right]_2-\left[\mathrm{PL}\right]\right)}=k_{\mathrm{on}}\left(\left[\mathrm{PL}\right]_1-\left[\mathrm{PL}\right]_2\right)t+\mathrm{C}_4

When :math:`t=0`, :math:`[\mathrm{PL}]=0`. Thus,

.. math::

    \ln{\left(\left[\mathrm{PL}\right]_1-0\right)}-\ln{\left(\left[\mathrm{PL}\right]_2-0\right)}=k_{\mathrm{on}}\left(\left[\mathrm{PL}\right]_1-\left[\mathrm{PL}\right]_2\right)\times0+\mathrm{C}_4

:math:`\mathrm{C}_4=\ln{\frac{\left[\mathrm{PL}\right]_1}{\left[\mathrm{PL}\right]_2}}`, substituting the above equation gives

.. math::

    \ln{\left(\left[\mathrm{PL}\right]_1-\left[\mathrm{PL}\right]\right)}-\ln{\left(\left[\mathrm{PL}\right]_2-\left[\mathrm{PL}\right]\right)}=k_{\mathrm{on}}\left(\left[\mathrm{PL}\right]_1-\left[\mathrm{PL}\right]_2\right)t+\ln{\frac{\left[\mathrm{PL}\right]_1}{\left[\mathrm{PL}\right]_2}}

.. math:: :label: ap_1.3

    t=\frac{1}{k_{\mathrm{on}}\left(\left[\mathrm{PL}\right]_1-\left[\mathrm{PL}\right]_2\right)}\ln{\frac{\left[\mathrm{PL}\right]_2\left(\left[\mathrm{PL}\right]_1-\left[\mathrm{PL}\right]\right)}{\left[\mathrm{PL}\right]_1\left(\left[\mathrm{PL}\right]_2-\left[\mathrm{PL}\right]\right)}}

Thermodynamics of the second-order binding process
--------------------------------------------------

The equilibrium state of the binding process is described by the following equation,

.. math::

    K_\mathrm{d}=\frac{\left[\mathrm{P}\right]_\mathrm{eq}\left[\mathrm{L}\right]_\mathrm{eq}}{\left[\mathrm{PL}\right]_\mathrm{eq}}

where the :math:`[\mathrm{P}]_\mathrm{eq}`, :math:`[\mathrm{L}]_\mathrm{eq}`, and :math:`[\mathrm{PL}]_\mathrm{eq}` is the equilibrium concentration of the unbound protein and ligand and
the equilibrium concentration of the protein-ligand complex. The total protein and ligand concentration,
:math:`[\mathrm{P}]_0` and :math:`[\mathrm{L}]_0`, and the equilibrium concentration of the unbound protein and ligand :math:`[\mathrm{P}]_\mathrm{eq}` and :math:`[\mathrm{L}]_\mathrm{eq}` are related by

.. math::

    \left[\mathrm{P}\right]_\mathrm{eq}=\left[\mathrm{P}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}

    \left[\mathrm{L}\right]_\mathrm{eq}=\left[\mathrm{L}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}

Thus,

.. math::

    K_\mathrm{d}=\frac{\left(\left[\mathrm{P}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}\right)\left(\left[\mathrm{L}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}\right)}{\left[\mathrm{PL}\right]_\mathrm{eq}}

    K_\mathrm{d}{\left[\mathrm{PL}\right]_\mathrm{eq}}=\left(\left[\mathrm{P}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}\right)\left(\left[\mathrm{L}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}\right)

.. math:: :label: ap_1.4

    {(\left[\mathrm{PL}\right]_\mathrm{eq})}^2-\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0-K_\mathrm{d}\right)\left[\mathrm{PL}\right]_\mathrm{eq}+\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_0=0


This equation is the same as the right-hand side of equation 2, its coefficients and discriminant are

.. math::

    a=1,\ b=-\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0+\frac{k_\mathrm{off}}{k_\mathrm{on}}\right)=-\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0+K_\mathrm{d}\right),\ c=\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_0

    \Delta=b^2-4ac=\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0+K_\mathrm{d}\right)^2-4\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_0=\left(\left[\mathrm{P}\right]_0-\left[L\right]_0\right)^2+{K_\mathrm{d}}^2+2\left[\mathrm{P}\right]_0K_\mathrm{d}+2\left[\mathrm{L}\right]_0K_\mathrm{d}>0

.. math:: :label: ap_1.5

    \left[\mathrm{PL}\right]_\mathrm{eq1}=\frac{-b-\sqrt{b^2-4ac}}{2a}=\frac{\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0+K_\mathrm{d}\right)-\sqrt{\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0+K_\mathrm{d}\right)^2-4\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_0}}{2}

    \left[\mathrm{PL}\right]_\mathrm{eq2}=\frac{-b+\sqrt{b^2-4ac}}{2a}=\frac{\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0+K_\mathrm{d}\right)+\sqrt{\left(\left[\mathrm{P}\right]_0+\left[\mathrm{L}\right]_0+K_\mathrm{d}\right)^2-4\left[\mathrm{P}\right]_0\left[\mathrm{L}\right]_0}}{2}

From the expressions of :math:`[\mathrm{PL}]_\mathrm{eq1}` and :math:`[\mathrm{PL}]_\mathrm{eq2}`, we know that :math:`0>[\mathrm{PL}]_\mathrm{eq1}>[\mathrm{PL}]_\mathrm{eq2}`. We assume that
the :math:`[\mathrm{PL}]` increases from 0 to the first real root to make the quadratic equation zero. So, :math:`[\mathrm{PL}]_\mathrm{eq1}`
should be the equilibrium concentration of the protein-ligand complex. It should be noted that the
thermodynamic process of the binding reaction is independent of the initial concentrations of the binding species.
Hence, we can calculate the total concentration of the protein and ligand from the initial binding system and use
the approaches presented here to calculate the equilibrium concentration of the protein, ligand, and protein-ligand complex.

Kinetics of the pseudo-first-order binding process
--------------------------------------------------

The binding of the protein and ligand is described by the following equation,

.. math::

    \mathrm{P} + \mathrm{L} \mathrel{\mathop{\rightleftarrows}^{k_{\mathrm{on}}}_{k_{\mathrm{off}}}} \mathrm{PL}

If :math:`\left[\mathrm{L}\right]_0\gg\left[\mathrm{P}\right]_0`, the binding process is effectively first-order since the
:math:`[\mathrm{L}]` is hardly affected by the binding process, then the equation can be transformed to

.. math::

    \mathrm{P}\mathrel{\mathop{\rightleftarrows}^{k_{\mathrm{on}}\times\mathrm{L}_0}_{k_{\mathrm{off}}}} \mathrm{PL}

The concentration of the protein-ligand complex changes during the binding process and

.. math::

    \frac{d\left[\mathrm{PL}\right]}{dt}=k_\mathrm{on}\left[\mathrm{P}\right]\left[\mathrm{L}\right]-k_\mathrm{off}\left[\mathrm{PL}\right]

During the binding process, :math:`\left[\mathrm{L}\right]\approx\left[\mathrm{L}\right]_0`, then

.. math::

    \frac{d\left[\mathrm{PL}\right]}{dt}=k_\mathrm{on}\left[\mathrm{P}\right]\left[\mathrm{L}\right]_0-k_\mathrm{off}\left[\mathrm{PL}\right]

The total protein concentration (:math:`[\mathrm{P}]_0`) and the concentration of the unbound protein (:math:`[\mathrm{P}]`) are related by

.. math::

    \left[\mathrm{P}\right]=\left[\mathrm{P}\right]_0-\left[\mathrm{PL}\right]

Thus,

.. math::

    \frac{d\left[\mathrm{PL}\right]}{dt}=k_\mathrm{on}\left[\mathrm{L}\right]_0\left(\left[\mathrm{P}\right]_0-\left[\mathrm{PL}\right]\right)-k_\mathrm{off}\left[\mathrm{PL}\right]

When :math:`[\mathrm{PL}]` increases from 0 to :math:`[\mathrm{PL}]_\mathrm{eq}`, :math:`\frac{d\left[\mathrm{PL}\right]}{dt}` decreases from the
maximum value to 0. :math:`\frac{d\left[\mathrm{PL}\right]}{dt}\rightarrow0` and :math:`\frac{d\left[\mathrm{PL}\right]}{dt}>0`,

.. math:: :label: ap_1.6

    \frac{d\left[\mathrm{PL}\right]}{dt}=k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)\left[\mathrm{PL}\right]>0

    \frac{d\left\{k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)\left[\mathrm{PL}\right]\right\}}{k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)\left[\mathrm{PL}\right]}=-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)dt

    \ln{\left|k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)\left[\mathrm{PL}\right]\right|}+\mathrm{C}_1=-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)t+\mathrm{C}_2

Because the right-hand side of equation 6 is greater than 0,

.. math::

    \ln{\left\{k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)\left[\mathrm{PL}\right]\right\}}=-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)t+\mathrm{C}_3

When :math:`t=0,\left[\mathrm{PL}\right]=0`. Thus, :math:`\mathrm{C}_3=\ln{\left(k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0\right)}`, substituting the above equation gives

.. math:: :label: ap_1.7

    \ln{\left\{k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)\left[\mathrm{PL}\right]\right\}}=-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)t+\ln{\left(k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0\right)}

    \left[\mathrm{PL}\right]=\frac{k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0}{k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}}\left(1-e^{-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)t}\right)

When :math:`t\rightarrow+\infty`,

.. math:: :label: ap_1.8

    \left[\mathrm{PL}\right]_{eq}=\frac{k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0}{k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}}=\frac{\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0}{\left[\mathrm{L}\right]_0+K_\mathrm{d}}

Substitute equation 8 with equation 7 gives

.. math:: :label: ap_1.9

    \left[\mathrm{PL}\right]=\left[\mathrm{PL}\right]_\mathrm{eq}\left(1-e^{-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)t}\right)

Define the observation rate constant :math:`k_\mathrm{obs}` by

.. math:: :label: ap_1.10

    k_\mathrm{obs}=k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}

Thus,

.. math:: :label: ap_1.11

    \left[\mathrm{PL}\right]=\left[\mathrm{PL}\right]_\mathrm{eq}\left(1-e^{-k_\mathrm{obs}t}\right)

Transform equation 7 to

.. math::

    \ln{\frac{k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)\left[\mathrm{PL}\right]}{k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0}}=-\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)t

.. math:: :label: ap_1.12

    t=\frac{-1}{k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}}\ln{\left(1-\frac{\left(k_\mathrm{on}\left[\mathrm{L}\right]_0+k_\mathrm{off}\right)\left[\mathrm{PL}\right]}{k_\mathrm{on}\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0}\right)}

Thermodynamics of the pseudo-first-order binding process
--------------------------------------------------------

The equilibrium state of the binding process is described by the following equation,

.. math::

    K_\mathrm{d}=\frac{\left[\mathrm{P}\right]_\mathrm{eq}\left[\mathrm{L}\right]_\mathrm{eq}}{\left[\mathrm{PL}\right]_\mathrm{eq}}

where :math:`[\mathrm{P}]_\mathrm{eq}`, :math:`[\mathrm{L}]_\mathrm{eq}`, and :math:`[\mathrm{PL}]_\mathrm{eq}` are the equilibrium concentrations of the unbound
protein and ligand and the equilibrium concentration of the protein-ligand complex. If
:math:`\left[\mathrm{L}\right]_0\gg\left[\mathrm{P}\right]_0`, or strictly speaking, :math:`\left[\mathrm{L}\right]_0\gg[\mathrm{PL}]_\mathrm{eq}`,
:math:`\left[\mathrm{L}\right]_\mathrm{eq}=\left[\mathrm{L}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}\approx\left[\mathrm{L}\right]_0`. Hence, the :math:`[\mathrm{L}]_0`
is hardly affected by the binding process, then the equation can be transformed to

.. math::

    K_\mathrm{d}=\frac{\left[\mathrm{P}\right]_\mathrm{eq}\left[\mathrm{L}\right]_0}{\left[\mathrm{PL}\right]_\mathrm{eq}}

The total protein concentration (:math:`[\mathrm{P}]_0`) and the equilibrium concentration of the unbound protein (:math:`[\mathrm{P}]_\mathrm{eq}`) are related by

.. math::

    \left[\mathrm{P}\right]_\mathrm{eq}=\left[\mathrm{P}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}

Thus,

.. math::

    K_\mathrm{d}=\frac{\left(\left[\mathrm{P}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}\right)\left[\mathrm{L}\right]_0}{\left[\mathrm{PL}\right]_\mathrm{eq}}

    K_\mathrm{d}{\left[\mathrm{PL}\right]_\mathrm{eq}}=\left(\left[\mathrm{P}\right]_0-\left[\mathrm{PL}\right]_\mathrm{eq}\right)\left[\mathrm{L}\right]_0

.. math:: :label: ap_1.13

    \left[\mathrm{PL}\right]_\mathrm{eq}=\frac{\left[\mathrm{L}\right]_0\left[\mathrm{P}\right]_0}{\left[\mathrm{L}\right]_0+K_\mathrm{d}}

In either the second-order binding process or the pseudo-first-order binding process, the :math:`[\mathrm{PL}]_\mathrm{eq}`
is the same calculated by the thermodynamic and kinetic approaches. But the equations calculated from the
second-order binding process are different from the equations calculated from the pseudo-first-order binding process.

Kinetics of the dissociation process
------------------------------------

The dissociation of the protein-ligand complex with no rebinding is described by the following equation,

.. math::

    \mathrm{PL}\overset{k_\mathrm{off}}{\rightarrow}\mathrm{P}+\mathrm{L}

The concentration of the protein-ligand complex changes during the dissociation process and

.. math::

    \frac{d\left[\mathrm{PL}\right]}{dt}=-k_\mathrm{off}\left[\mathrm{PL}\right]

    \frac{d\left[\mathrm{PL}\right]}{\left[\mathrm{PL}\right]}=-k_\mathrm{off}dt

    \ln{\left|\left[\mathrm{PL}\right]\right|}=-k_\mathrm{off}t+\mathrm{C}

    \ln{\left[\mathrm{PL}\right]}=-k_\mathrm{off}t+\mathrm{C}

When :math:`t=0, \left[\mathrm{PL}\right]=\left[\mathrm{PL}\right]_0`. Thus, :math:`\mathrm{C}=\ln{\left[\mathrm{PL}\right]_0}`, substituting the above equation gives

.. math::

    \ln{\left[\mathrm{PL}\right]}=-k_\mathrm{off}t+\ln{\left[\mathrm{PL}\right]_0}


.. math:: :label: ap_1.14

    t=-\frac{1}{k_\mathrm{off}}\ln{\frac{\left[\mathrm{PL}\right]}{\left[\mathrm{PL}\right]_0}}
