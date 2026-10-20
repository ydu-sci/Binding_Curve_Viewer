==========
Appendix 3
==========

Calculation of IC50 and Apparent IC50
=====================================

At equilibrium, the concentrations of all chemical species in the competitive
binding system are governed by Equations :eq:`ap3_2`–:eq:`ap3_6`, where
:math:`[\mathrm{P}]_0`, :math:`[\mathrm{L}]_0`, and :math:`[\mathrm{I}]_0`
represent the total concentrations of protein, ligand, and inhibitor,
respectively. Assuming that :math:`[\mathrm{P}]_0`, :math:`[\mathrm{L}]_0`,
:math:`K_\mathrm{d}`, and :math:`K_\mathrm{i}` are known, Equations
:eq:`ap3_7` and :eq:`ap3_8` can be derived from Equations :eq:`ap3_2` and
:eq:`ap3_4`, and from Equations :eq:`ap3_3` and :eq:`ap3_5`, respectively.

When the protein–ligand complex is inhibited by 50%, the equilibrium complex
concentration is defined as
:math:`[\mathrm{PL}]_{\mathrm{eq\text{-}50}} = 0.5 \times [\mathrm{PL}]_{\mathrm{eq\text{-}0}}`,
where :math:`[\mathrm{PL}]_{\mathrm{eq\text{-}0}}` denotes the equilibrium
concentration of the protein–ligand complex in the absence of inhibitor.
:math:`[\mathrm{PL}]_{\mathrm{eq\text{-}0}}` can be calculated using
Equation :eq:`ap3_1`.

Given :math:`[\mathrm{PL}]_{\mathrm{eq\text{-}50}}`, the equilibrium
concentration of free protein :math:`[\mathrm{P}]_{\mathrm{eq\text{-}50}}`
can be determined by Equation :eq:`ap3_7`, and the equilibrium concentration
of the protein–inhibitor complex
:math:`[\mathrm{PI}]_{\mathrm{eq\text{-}50}}` can be determined by Equation
:eq:`ap3_6`. The apparent IC50 is the :math:`[\mathrm{I}]_0` that elicits a
50% decrease in :math:`[\mathrm{PL}]_{\mathrm{eq\text{-}0}}`, and it can be
calculated from Equation :eq:`ap3_8`. Consequently, the IC50, defined as the
free inhibitor concentration :math:`[\mathrm{I}]` when the protein–ligand
complex is displaced by 50%, can be calculated using Equation :eq:`ap3_5`.

By substituting Equation :eq:`ap3_8` into Equation :eq:`ap3_5`, we obtain
Equation :eq:`ap3_9`. Since :math:`[\mathrm{PI}]_{\mathrm{eq\text{-}50}}`
and :math:`[\mathrm{P}]_{\mathrm{eq\text{-}50}}` are fully determined by
:math:`[\mathrm{PL}]_{\mathrm{eq\text{-}50}}`, :math:`[\mathrm{P}]_0`,
:math:`[\mathrm{L}]_0`, and :math:`K_\mathrm{d}`, the IC50 (i.e.,
:math:`[\mathrm{I}]`) is proportional to :math:`K_\mathrm{i}`. However,
the apparent IC50 is not directly proportional to :math:`K_\mathrm{i}`.

.. math::
   :label: ap3_1

   [\mathrm{PL}]_\mathrm{eq} =
   \frac{\left([\mathrm{P}]_0 + [\mathrm{L}]_0 + K_\mathrm{d}\right)
   - \sqrt{\left([\mathrm{P}]_0 + [\mathrm{L}]_0 + K_\mathrm{d}\right)^2
   - 4[\mathrm{P}]_0[\mathrm{L}]_0}}{2}

.. math::
   :label: ap3_2

   K_\mathrm{d} = \frac{[\mathrm{P}][\mathrm{L}]}{[\mathrm{PL}]}

.. math::
   :label: ap3_3

   K_\mathrm{i} = \frac{[\mathrm{P}][\mathrm{I}]}{[\mathrm{PI}]}

.. math::
   :label: ap3_4

   [\mathrm{L}] = [\mathrm{L}]_0 - [\mathrm{PL}]

.. math::
   :label: ap3_5

   [\mathrm{I}] = [\mathrm{I}]_0 - [\mathrm{PI}]

.. math::
   :label: ap3_6

   [\mathrm{PI}] = [\mathrm{P}]_0 - [\mathrm{P}] - [\mathrm{PL}]

.. math::
   :label: ap3_7

   [\mathrm{P}] = \frac{K_\mathrm{d}\,[\mathrm{PL}]}{[\mathrm{L}]_0 - [\mathrm{PL}]}

.. math::
   :label: ap3_8

   [\mathrm{I}]_0 = \frac{[\mathrm{PI}]\left(K_\mathrm{i} + [\mathrm{P}]\right)}{[\mathrm{P}]}

.. math::
   :label: ap3_9

   [\mathrm{I}] = \frac{[\mathrm{PI}]\,K_\mathrm{i}}{[\mathrm{P}]}
