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
from bokeh.plotting import figure, output_file, show
from bokeh.models import Spinner, TabPanel, Tabs, CustomJS, ColumnDataSource, TapTool, Circle, \
    LabelSet, DataTable, TableColumn, CDSView, BooleanFilter, HTMLTemplateFormatter, Div, Scatter
from bokeh.models.tools import HoverTool, CrosshairTool
from bokeh.layouts import column, row, layout
from bokeh.models import ColorBar, LinearColorMapper, Rect
from bokeh.models import BasicTicker, PrintfTickFormatter
from bokeh.transform import linear_cmap



def dissociation_constant(pt_const=100.0, kd_const=100.0, plot_points=1000, shown_points_num=971):
    r"""The function ``dissociation_constant`` generates the saturation curves of protein-ligand
    complex. Here, we specify that the protein is the component in fixed concentration. The
    ligand can be any molecules, e.g., DNA, RNA, or peptides, that binds to the target protein. The
    dissociation constant, :math:`K_\mathrm{d}`, is defined by the equation

    .. math:: :label: eq_1.1

        K_\mathrm{d}=\frac{[\mathrm{P}][\mathrm{L}]}{[\mathrm{PL}]}

    where :math:`[\mathrm{P}]` is the concentration of the free protein, :math:`[\mathrm{L}]` is
    the concentration of the free ligand, and :math:`[\mathrm{PL}]` is the concentration of the
    protein-ligand complex.

    According to the parameters of :math:`K_\mathrm{d}` and concentration of the total protein
    :math:`[\mathrm{P}]_{\mathrm{total}}`, the following equation derived from equation 1 is used
    to generate the saturation curve by a step size sequence of :math:`[\mathrm{PL}]`.

    .. math:: :label: eq_1.2

        [\mathrm{L}]_{\mathrm{total}}=\frac{([\mathrm{P}]_{\mathrm{total}}+K_\mathrm{d}-[\mathrm{PL}])[\mathrm{PL}]}{[\mathrm{P}]_{\mathrm{total}}-[\mathrm{PL}]}

    Parameters
    ----------
    pt_const : float
        The concentration of the total protein (:math:`[\mathrm{P}]_{\mathrm{total}}`) in the unit of nM (default value is 100).
    kd_const : float
        The value of the :math:`K_\mathrm{d}` in the unit of nM (default value is 100).
    plot_points : int
        The number of points (i.e., smoothness) of lines (default value is 1000).
    shown_points_num : int
        The number of points in lines shown in the figure (default value is 971). To avoid the
        long tail before proteins are fully occupied by high concentration of ligand, by default,
        binding curves end when 97.0% of total protein are bound. ``shown_points_num`` should be
        less than ``plot_points``.

    Returns
    -------
    Tuple
        A Tuple of widgets and figures, including ``p1_header``, ``config_col``, ``tab_conc``,
        ``tab_frac`` can be used to generate the web application.
    """
    # pt = [P]total, lt = [L]total, pl = [PL], pf = [P]free, lf = [L]free
    # kd = pf*lf / pl => kd = (pt-pl) * (lt-pl) / pl
    # => lt = (pt+kd-pl) * pl / (pt-pl)
    num_array = np.arange(0, plot_points, 1) / plot_points
    kd = np.full(shape=shown_points_num, fill_value=kd_const)
    pl_tmp = pt_const * num_array
    # define the range of [PL] (Y-axis)
    pl = pl_tmp[: shown_points_num]
    pt = np.full(shape=shown_points_num, fill_value=pt_const)
    lt = (pt + kd - pl) * pl / (pt - pl)
    lf = lt - pl
    lf[0] = 0  # change float 0.0000 to int 0
    lt[0] = 0  # change float 0.0000 to int 0
    app_kd = lt[500]
    lt_95 = lt[950]
    lf_95 = lf[950]
    fraction_bp = pl / pt  # fraction of bound protein
    # Figure generation - P1-1
    source = ColumnDataSource(data=dict(pt=pt, pl=pl, kd=kd, lt=lt, lf=lf, fraction_bp=fraction_bp))
    plot_linear = figure(width=400, height=400, tools=[], output_backend='svg')
    plot_linear.toolbar.logo = None
    plot_linear.xaxis.axis_label = "Concentration of Free/Total Ligand (nM)"
    plot_linear.yaxis.axis_label = "Concentration of Bound Protein ([PL], nM)"
    plot_linear.xaxis.axis_label_text_font_style = 'normal'
    plot_linear.yaxis.axis_label_text_font_style = 'normal'
    plot_linear.axis.minor_tick_in = -3
    plot_linear.axis.minor_tick_out = 6
    lf_line = plot_linear.line('lf', 'pl', source=source, legend_label='[L]free', line_width=2, color='#66c2a5')
    lt_line = plot_linear.line('lt', 'pl', source=source, legend_label='[L]total', line_width=2, color='#fc8d62',
                               name='lt_line')
    plot_linear.legend.location = 'bottom_right'
    htools_lt_line = HoverTool(tooltips=[
        ('[L]free (nM)', '@lf{%0.2e}'),
        ('[L]total (nM)', '@lt{%0.2e}'),
        ('[PL] (nM)', '@pl{%0.2e}'), ],
        formatters={'@lf': 'printf', '@lt': 'printf', '@pl': 'printf'},
        mode='hline')
    htools_lt_line.renderers = [lt_line]
    ctools_lt_line = CrosshairTool(dimensions='width', line_alpha=0.5, line_color='grey', line_width=1)
    plot_linear.add_tools(htools_lt_line, ctools_lt_line)
    tab_conc1 = TabPanel(child=plot_linear, title='Linear X-axis')
    # Figure generation - P2-1
    plot_linear2 = figure(width=400, height=400, tools=[], output_backend='svg')
    plot_linear2.toolbar.logo = None
    plot_linear2.xaxis.axis_label = "Concentration of Free/Total Ligand (nM)"
    plot_linear2.yaxis.axis_label = r"Fraction of Bound Protein ($$\text{[PL]}/\text{[P]}_0$$)"
    plot_linear2.xaxis.axis_label_text_font_style = 'normal'
    plot_linear2.yaxis.axis_label_text_font_style = 'normal'
    plot_linear2.axis.minor_tick_in = -3
    plot_linear2.axis.minor_tick_out = 6
    lf_line2 = plot_linear2.line('lf', 'fraction_bp', source=source, legend_label='[L]free', line_width=2,
                                 color='#66c2a5')
    lt_line2 = plot_linear2.line('lt', 'fraction_bp', source=source, legend_label='[L]total', line_width=2,
                                 color='#fc8d62', name='lt_line2')
    plot_linear2.legend.location = 'bottom_right'
    htools_lt_line2 = HoverTool(tooltips=[
        ('[L]free (nM)', '@lf{%0.2e}'),
        ('[L]total (nM)', '@lt{%0.2e}'),
        ('[PL]/[P]total', '@fraction_bp{(0.00)}'), ],
        formatters={'@lf': 'printf', '@lt': 'printf'},
        mode='hline')
    htools_lt_line2.renderers = [lt_line2]
    ctools_lt_line2 = CrosshairTool(dimensions='width', line_alpha=0.5, line_color='grey', line_width=1)
    plot_linear2.add_tools(htools_lt_line2, ctools_lt_line2)
    tab_frac1 = TabPanel(child=plot_linear2, title='Linear X-axis')
    # Figure generation - P1-2
    plot_log = figure(width=400, height=400, x_axis_type='log', tools=[], output_backend='svg')
    plot_log.toolbar.logo = None
    plot_log.xaxis.axis_label = "Concentration of Free/Total Ligand (nM)"
    plot_log.yaxis.axis_label = "Concentration of Bound Protein ([PL], nM)"
    plot_log.xaxis.axis_label_text_font_style = 'normal'
    plot_log.yaxis.axis_label_text_font_style = 'normal'
    plot_log.axis.minor_tick_in = -3
    plot_log.axis.minor_tick_out = 6
    lf_line_log = plot_log.line('lf', 'pl', source=source, legend_label='[L]free', line_width=2, color='#66c2a5')
    lt_line_log = plot_log.line('lt', 'pl', source=source, legend_label='[L]total', line_width=2, color='#fc8d62',
                                name='lt_line_log')
    plot_log.legend.location = 'bottom_right'
    htools_lt_line_log = HoverTool(tooltips=[
        ('[L]free (nM)', '@lf{%0.2e}'),
        ('[L]total (nM)', '@lt{%0.2e}'),
        ('[PL] (nM)', '@pl{%0.2e}'), ],
        formatters={'@lf': 'printf', '@lt': 'printf', '@pl': 'printf'},
        mode='hline')
    htools_lt_line_log.renderers = [lt_line_log]
    ctools_lt_line_log = CrosshairTool(dimensions='width', line_alpha=0.5, line_color='grey', line_width=1)
    plot_log.add_tools(htools_lt_line_log, ctools_lt_line_log)
    tab_conc2 = TabPanel(child=plot_log, title='Logarithmic X-axis')
    # Figure generation - P2-2
    plot_log2 = figure(width=400, height=400, x_axis_type='log', tools=[], output_backend='svg')
    plot_log2.toolbar.logo = None
    plot_log2.xaxis.axis_label = "Concentration of Free/Total Ligand (nM)"
    plot_log2.yaxis.axis_label = r"Fraction of Bound Protein ($$\text{[PL]}/\text{[P]}_0$$)"
    plot_log2.xaxis.axis_label_text_font_style = 'normal'
    plot_log2.yaxis.axis_label_text_font_style = 'normal'
    plot_log2.axis.minor_tick_in = -3
    plot_log2.axis.minor_tick_out = 6
    lf_line_log2 = plot_log2.line('lf', 'fraction_bp', source=source, legend_label='[L]free', line_width=2,
                                  color='#66c2a5')
    lt_line_log2 = plot_log2.line('lt', 'fraction_bp', source=source, legend_label='[L]total', line_width=2,
                                  color='#fc8d62', name='lt_line_log2')
    plot_log2.legend.location = 'bottom_right'
    htools_lt_line_log2 = HoverTool(tooltips=[
        ('[L]free (nM)', '@lf{%0.2e}'),
        ('[L]total (nM)', '@lt{%0.2e}'),
        ('[PL]/[P]total', '@fraction_bp{(0.00)}'), ],
        formatters={'@lf': 'printf', '@lt': 'printf'},
        mode='hline')
    htools_lt_line_log2.renderers = [lt_line_log2]
    ctools_lt_line_log2 = CrosshairTool(dimensions='width', line_alpha=0.5, line_color='grey', line_width=1)
    plot_log2.add_tools(htools_lt_line_log2, ctools_lt_line_log2)
    tab_frac2 = TabPanel(child=plot_log2, title='Logarithmic X-axis')
    # Callbacks and widgets
    p1_header = Div(text='<h1 style="display:inline; margin-top: 0">Determination of Dissociation Constant '
                         +'(<i>K</i><sub>d</sub>)</h1>'
                    + '<span>&nbsp;<a href="/binding-curve-viewer/">[Back to Home]</a></span>',
                    height=40, width=800, height_policy='auto', margin=(5, 5, 5, 5), )
    data_table = Div(text = 'Apparent <b><i>K</i><sub>d</sub></b> = '+ '{:.5f}'.format(app_kd) +' nM; <b><i>K</i><sub>d</sub></b> = '+ '{:.5f}'.format(kd_const) +' nM;<br>'
             + 'Apparent <b><i>K</i><sub>d</sub></b> / <b><i>K</i><sub>d</sub></b> = '+ '{:.2f}'.format(app_kd / kd_const) +';<br>'
             + '<b>[L]<sub>free</sub></b> = '+ '{:.5f}'.format(lf_95) +' nM when 95% Protein is bound;<br>'
             + '<b>[L]<sub>total</sub></b> = '+ '{:.5f}'.format(lt_95) +' nM when 95% Protein is bound.',
        # This width is used to control the first column.
        width=300, margin=(15, 5, 5, 5))
    callback_kd = CustomJS(args=dict(source=source, shown_points_num=shown_points_num, data_table=data_table),
                           code="""
        const data = source.data;
        const user_kd = cb_obj.value;
        const pt = data['pt'];
        const pl = data['pl'];
        const kd = data['kd'];
        const lt = data['lt'];
        const lf = data['lf'];
        const fraction_bp = data['fraction_bp'];
        var app_kd_html;
        var kd_html;
        var lt_html;
        var lf_html;
        for (var i = 0; i < shown_points_num; i++) {
            kd[i] = user_kd;
            lt[i] = (pt[i]+user_kd-pl[i]) * pl[i] / (pt[i]-pl[i]);
            lf[i] = lt[i] - pl[i];
            if (i==500) {
                kd_html = lf[i];
                app_kd_html = lt[i];
            }
            if (i==950) {
                lt_html = lt[i];
                lf_html = lf[i];
            }
        }
        source.change.emit();
        const ratio_app_kd = app_kd_html / kd_html;
        data_table.text='Apparent <b><i>K</i><sub>d</sub></b> = '+app_kd_html.toFixed(5)
                        +' nM; <b><i>K</i><sub>d</sub></b> = '+kd_html.toFixed(5)+' nM;<br>' 
                        +'Apparent <b><i>K</i><sub>d</sub></b> / <b><i>K</i><sub>d</sub></b> = '
                        +ratio_app_kd.toFixed(2)+';<br>' + '<b>[L]<sub>free</sub></b> = '
                        +lf_html.toFixed(5)+' nM when 95% Protein is bound;<br>' 
                        +'<b>[L]<sub>total</sub></b> = '+lt_html.toFixed(5)+' nM when 95% Protein is bound.';
    """)
    callback_pt = CustomJS(
        args=dict(source=source, shown_points_num=shown_points_num, num_array=num_array, data_table=data_table),
        code="""
        const data = source.data;
        const user_pt = cb_obj.value;
        const pt = data['pt'];
        const pl = data['pl'];
        const kd = data['kd'];
        const lt = data['lt'];
        const lf = data['lf'];
        const fraction_bp = data['fraction_bp']
        var app_kd_html;
        var kd_html;
        var lt_html;
        var lf_html;
        for (var i=0; i<shown_points_num; i++) {
            pl[i] = user_pt * num_array[i];
            pt[i] = user_pt;
            lt[i] = (user_pt+kd[i]-pl[i]) * pl[i] / (user_pt-pl[i]);
            lf[i] = lt[i] - pl[i];
            fraction_bp[i] = pl[i] / user_pt;
            if (i==500) {
                kd_html = lf[i];
                app_kd_html = lt[i];
            }
            if (i==950) {
                lt_html = lt[i];
                lf_html = lf[i];
            }
        };
        source.change.emit();
        const ratio_app_kd = app_kd_html / kd_html;
        data_table.text='Apparent <b><i>K</i><sub>d</sub></b> = '+app_kd_html.toFixed(5)
                        +' nM; <b><i>K</i><sub>d</sub></b> = '+kd_html.toFixed(5)+' nM;<br>' 
                        +'Apparent <b><i>K</i><sub>d</sub></b> / <b><i>K</i><sub>d</sub></b> = '
                        +ratio_app_kd.toFixed(2)+';<br>' +'<b>[L]<sub>free</sub></b> = '
                        +lf_html.toFixed(5)+' nM when 95% Protein is bound;<br>' 
                        +'<b>[L]<sub>total</sub></b> = '+lt_html.toFixed(5)+' nM when 95% Protein is bound.';
    """)

    div_kd = Div(text="<p>Equilibrium Dissociation Constant (<b><i>K</i><sub>d</sub></b>, nM):" +
                      "<br>0.001 ~ 50000 nM</p>", width=300, height=35, margin=(0, 0, 0, 5))
    spinner_kd = Spinner(low=0.001, high=50000, value=kd_const, step=1, title='', width=240, name='spinner_kd')
    spinner_kd.js_on_change('value', callback_kd)
    div_pt = Div(text="<p>Concentration of Total Fixed Component (e.g., Protein, <b>[P]<sub>total</sub></b>, nM):" +
                      "<br>0.001 ~ 50000 nM</p>", width=300, height=48, margin=(5, 0, 0, 5))
    spinner_pt = Spinner(low=0.001, high=50000, value=pt_const, step=1, title='', width=240, name='spinner_pt')
    spinner_pt.js_on_change('value', callback_pt)
    config_col = column(div_kd, spinner_kd, div_pt, spinner_pt, data_table)
    tab_conc = Tabs(tabs=[tab_conc1, tab_conc2])
    tab_frac = Tabs(tabs=[tab_frac1, tab_frac2])
    return (p1_header, config_col, tab_conc, tab_frac)


def kinetics_association_dissociation(p0=100, l0=150, koff=1e-2, kon=1e5, plot_points=1000, shown_points_num=991):
    r"""The function ``kinetics_association_dissociation`` generates the association and dissociation
    curves of protein-ligand complex. The association kinetics is described by the pseudo-first
    order process and second-order process. The ligand can be any molecules, e.g., DNA, RNA, or
    peptides, that binds to the target protein. The dissociation constant, :math:`K_\mathrm{d}`,
    is defined by the equation

    .. math:: :label: eq_2.1

        K_\mathrm{d}=\frac{k_{\mathrm{off}}}{k_{\mathrm{on}}}

    where :math:`k_{\mathrm{off}}` is the dissociation rate constant, :math:`k_{\mathrm{on}}` is the association
    rate constant.

    In the second-order binding process, the relationship between the concentration of the protein-ligand
    complex and the time is

    .. math:: :label: eq_2.2

        t=\frac{1}{k_{\mathrm{on}}([\mathrm{PL}]_1-[\mathrm{PL}]_2)}\ln\frac{[\mathrm{PL}]_2([\mathrm{PL}]_1-[\mathrm{PL}])}{[\mathrm{PL}]_1([\mathrm{PL}]_2-[\mathrm{PL}])}

    where the :math:`[\mathrm{PL}]_1` and :math:`[\mathrm{PL}]_2` are the roots of a quadratic equation and they
    can be calculated from the :math:`[\mathrm{P}]_0`, :math:`[\mathrm{L}]_0`, and :math:`K_\mathrm{d}` (:math:`[\mathrm{P}]_0`
    and :math:`[\mathrm{L}]_0` is the initial concentration of the protein and ligand, see :ref:`Kinetics of the
    second-order binding process<kinetics_2nd>` in Appendix 1 for their expressions).

    In the pseudo-first-order binding process, the relationship between the concentration of the
    protein-ligand complex and the time is

    .. math:: :label: eq_2.3

        [\mathrm{PL}]=\frac{k_{\mathrm{on}}[\mathrm{L}]_0[\mathrm{P}]_0}{k_{\mathrm{on}}[\mathrm{L}]_0+k_{\mathrm{off}}}(1-\exp(-(k_{\mathrm{on}}[\mathrm{L}]_0+k_{\mathrm{off}})t))


    Parameters
    ----------
    p0 : float
        The initial concentration of the free protein in the unit of nM (default value is 100).
    l0 : float
        The initial concentration of the free ligand in the unit of nM (default value is 150).
    koff : float
        The value of the dissociation rate constant (:math:`k_{\mathrm{off}}`) in the unit of s\ :sup:`-1`
        (default value is 1e-2).
    kon : float
        The value of the association rate constant (:math:`k_{\mathrm{on}}`) in the unit of M\ :sup:`-1` s\ :sup:`-1`
        (default value is 1e5).
    plot_points : int
        The number of points (i.e., smoothness) of lines (default value is 1000).
    shown_points_num : int
        The number of points in lines shown in the figure (default value is 991). To avoid the
        long tail before equilibrium, the binding curve ends when 99.0% equilibrium is reached,
        shown_points_num should be less than plot_points.

    Returns
    -------
    Tuple
        A tuple of widgets and figures, including ``kin_header``, ``kin_config_col``, ``plot_kinetics_on``,
        ``plot_kinetics_disso`` can be used to generate the web application.
    """
    kin_data = dict(
        kon=[], koff=[], p0=[], l0=[],
        pl_asso_list=[], pseudo_pl_asso_list=[], time_list=[],
        pl_asso_frac_list_vs_pl_1=[], pl_asso_frac_list_vs_p0=[],
        pl_disso_list=[], pl_disso_time_list=[])
    # Pic 1's two horizontal annotation lines
    hline_data = dict(hline_x=[0, 0], hline_pl_1=[], hline_pl_eq_pseudo=[])
    # Pic 1's two annotation texts
    anno_text_data = dict(text_x=[0], text_y=[], text_y_pseudo=[], asso_text=[], asso_text_pseudo=[], disso_text=[])
    div_text_data = dict(t99=[], pseudo_t99=[], pt99_ratio=[], disso_t50=[], disso_t99=[])
    kin_data['p0'].extend([p0] * shown_points_num)
    kin_data['l0'].extend([l0] * shown_points_num)
    # In this function, the unit of kon is nM^-1 s^-1, converted from parameter kon's M-1 s-1
    # For example, the kinetics parameters of Kd 1 nM
    # kofflist = np.array([1e-1, 1e-3, 1e-6])       # s^-1
    # konlist = np.array([1e8, 1e6, 1e3]) * 1e-9    # convert M^-1 s^-1 to nM^-1 s^-1
    # konlist = np.array([1e-1, 1e-3, 1e-6])        # nM^-1 s^-1
    kon = kon * 1e-9  # convert the unit from M^-1 s^-1 to nM^-1 s^-1
    kd = koff / kon
    kin_data['kon'].extend([kon] * shown_points_num)
    kin_data['koff'].extend([koff] * shown_points_num)
    # a, b, c are the coefficients of quadratic equation
    # For more details, see the supplementary information of binding curve viewer paper
    a = 1
    b = 0 - (p0 + l0 + koff/kon)
    c = p0 * l0
    pl_1 = (0 - b - np.sqrt(b * b - 4 * a * c)) / (2 * a)  # pl_1 is [PL]eq, pl_1 < pl_2
    pl_2 = (0 - b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
    A = 1 / (pl_1 - pl_2)
    num_array = np.arange(0, plot_points, 1) / plot_points
    pl_asso_list = pl_1 * num_array[:shown_points_num]
    time_list = A / kon * np.log((pl_2 * (pl_1 - pl_asso_list)) / (pl_1 * (pl_2 - pl_asso_list)))
    pl_99perc = 0.99 * pl_1  # [PL] at 99% equilibrium
    t99 = A / kon * np.log((pl_2 * (pl_1 - pl_99perc)) / (pl_1 * (pl_2 - pl_99perc)))  # time for 99% equilibrium
    pl_asso_frac_list_vs_pl_1 = pl_asso_list / pl_1
    pl_asso_frac_list_vs_p0 = pl_asso_list / p0

    # Pseudo rate section l0>>p0
    # Kinetic and thermodynamic approaches generate the same Peq result, see paper for more details.
    pl_eq_pseudo = kon * l0 * p0 / (kon * l0 + koff)
    kobs = kon * l0 + koff
    # We used the time_list for pseudo binding curve to keep the both curves in the same time frame.
    # The following equations are for your information.
    # pseudo_pl_asso_list = p_eq_pseudo * num_array[:shown_points_num]
    # pseudo_time_list = ((-1)/kobs)*np.log(1-kobs*pseudo_pl_asso_list/(kon*p0*l0))
    pseudo_pl_asso_list = (kon * p0 * l0 / kobs) * (1 - np.exp(0 - kobs * time_list))
    pseudo_t99 = np.log(100) / kobs

    # Dissociation section
    disso_num_array = np.arange(plot_points, 0, -1) / plot_points
    pl_disso_list = pl_1 * disso_num_array[:shown_points_num]
    pl_disso_time_list = (-1 / koff * np.log(pl_disso_list / pl_1))
    disso_t99 = np.log(100) / koff
    disso_t50 = np.log(2) / koff
    pt99_ratio = pseudo_t99 / t99


    kin_data['time_list'] = time_list[:]
    kin_data['pl_asso_list'] = pl_asso_list[:]
    kin_data['pseudo_pl_asso_list'] = pseudo_pl_asso_list[:]
    kin_data['pl_asso_frac_list_vs_pl_1'] = pl_asso_frac_list_vs_pl_1[:]
    kin_data['pl_asso_frac_list_vs_p0'] = pl_asso_frac_list_vs_p0[:]
    hline_data['hline_x'][1] = t99 * 1.2
    kin_data['pl_disso_list'] = pl_disso_list[:]
    kin_data['pl_disso_time_list'] = pl_disso_time_list[:]
    hline_data['hline_pl_1'] = [pl_1, pl_1]
    hline_data['hline_pl_eq_pseudo'] = [pl_eq_pseudo, pl_eq_pseudo]
    anno_text_data['text_y'] = [pl_1]
    anno_text_data['text_y_pseudo'] = [pl_eq_pseudo]
    anno_text_data['asso_text'] = ['[PL]eq in 2nd-order process: %.2e nM' % (pl_1)]
    anno_text_data['asso_text_pseudo'] = ['[PL]eq in pseudo-1st-order process: %.2e nM' % (pl_eq_pseudo)]
    anno_text_data['disso_text'] = ['[PL]0 = %.2e nM' % (pl_1)]
    div_text_data['t99'] = [t99]
    div_text_data['pseudo_t99'] = [pseudo_t99]
    div_text_data['pt99_ratio'] = [pseudo_t99 / t99]
    div_text_data['disso_t50'] = [disso_t50]
    div_text_data['disso_t99'] = [disso_t99]
    kin_source = ColumnDataSource(data=kin_data)
    hline_source = ColumnDataSource(data=hline_data)
    anno_text_source = ColumnDataSource(data=anno_text_data)
    div_text_source = ColumnDataSource(data=div_text_data)

    # Figure generation of equilibria
    plot_kinetics_on = figure(width=580, height=600, y_range=(-0.05 * p0, 1.25 * p0), tools=[], margin=(0, 0, 0, 5),
                              output_backend='svg')
    plot_kinetics_on.toolbar.logo = None
    plot_kinetics_on.xaxis.axis_label = "Time (s)"
    plot_kinetics_on.yaxis.axis_label = "Concentration of Bound Protein ([PL], nM)"
    plot_kinetics_on.xaxis.axis_label_text_font_style = 'normal'
    plot_kinetics_on.yaxis.axis_label_text_font_style = 'normal'
    plot_kinetics_on.axis.minor_tick_in = -3
    plot_kinetics_on.axis.minor_tick_out = 6
    pseudo_on_line = plot_kinetics_on.line('time_list', 'pseudo_pl_asso_list', source=kin_source, color='#66c2a5',
                                           name='pseudo_on_line', legend_label='Pseudo-1st-order process',
                                           line_width=2)
    on_line = plot_kinetics_on.line('time_list', 'pl_asso_list', source=kin_source, color='#fc8d62', name='on_line',
                                    legend_label='2nd-order process', line_width=2, line_alpha=0.5)
    dash_line1 = plot_kinetics_on.line('hline_x', 'hline_pl_1', source=hline_source, line_alpha=0.8,
                                       line_dash='dotdash', color='#fc8d62')
    dash_line2 = plot_kinetics_on.line('hline_x', 'hline_pl_eq_pseudo', source=hline_source, line_alpha=0.8,
                                       line_dash='dotdash', color='#66c2a5', )
    # Annotation texts and lines.
    l1 = LabelSet(x='text_x', y='text_y', x_units='screen', text='asso_text', source=anno_text_source,
                  border_line_color=None, x_offset=0, y_offset=-20, border_line_alpha=1.0,
                  background_fill_color=None, background_fill_alpha=1.0)
    l2 = LabelSet(x='text_x', y='text_y_pseudo', x_units='screen', text='asso_text_pseudo', source=anno_text_source,
                  border_line_color=None, x_offset=0, y_offset=5, border_line_alpha=1.0,
                  background_fill_color=None, background_fill_alpha=1.0)
    plot_kinetics_on.add_layout(l1)
    plot_kinetics_on.add_layout(l2)
    kin_vtools = HoverTool(tooltips=[
        ('Tims (s)', '@time_list{(0.00)}'),
        ('[PL] in pseudo-1st-order process (nM)', '@pseudo_pl_asso_list{%0.1e}'),
        ('[PL] in 2nd-order process (nM)', '@pl_asso_list{%0.1e}'),
        ('[PL]/[PL]eq in 2nd-order process', '@pl_asso_frac_list_vs_pl_1{%0.1e}'),
        ('[PL]/[P]0 in 2nd-order process', '@pl_asso_frac_list_vs_p0{%0.1e}')],
        formatters={'@pseudo_pl_asso_list': 'printf', '@pl_asso_list': 'printf'},
        mode='vline')
    kin_vtools.renderers = [on_line]
    kin_ctools = CrosshairTool(dimensions='height', line_alpha=0.5, line_color='grey', line_width=1)
    plot_kinetics_on.add_tools(kin_vtools, kin_ctools)
    plot_kinetics_on.legend.location = 'top_left'

    kin_data_div = Div(
        text='Time to Reach 99% Equilibrium in Pseudo-1st-order Process '
             + '(<b><i>t</i><sub>pseudo-0.99</sub></b>) = ' + '%.2e' % div_text_data['pseudo_t99'][0]
             + ' s;<br>' + 'Time to Reach 99% Equilibrium in 2nd-order Process '
             + '(<b><i>t</i><sub>0.99</sub></b>) = ' + '%.2e' % div_text_data['t99'][0]
             + ' s;<br>' + '<b><i>t</i><sub>pseudo-0.99</sub></b>/<b><i>t</i><sub>0.99</sub></b>'
             + ' = %.2f' % div_text_data['pt99_ratio'][0] + ';<br>'
             + 'Half-life of Dissociation (ln(2)/koff) = %.2e' % div_text_data['disso_t50'][0]
             + ' s;<br>' + 'Time to Reach 99% Dissociation = {:.2e}'.format(div_text_data['disso_t99'][0]) + ' s.',
        width=300, margin=(0, 0, 5, 5))

    # Figure generation of dissociation
    plot_kinetics_disso = figure(width=580, height=600, tools=[], output_backend='svg',
                                 y_range=(-0.05 * kin_data['pl_disso_list'][0], 1.1 * kin_data['pl_disso_list'][0]))
    plot_kinetics_disso.toolbar.logo = None
    plot_kinetics_disso.xaxis.axis_label = "Time (s)"
    plot_kinetics_disso.yaxis.axis_label = "Concentration of Bound Protein ([PL], nM)"
    plot_kinetics_disso.xaxis.axis_label_text_font_style = 'normal'
    plot_kinetics_disso.yaxis.axis_label_text_font_style = 'normal'
    plot_kinetics_disso.axis.minor_tick_in = -3
    plot_kinetics_disso.axis.minor_tick_out = 6
    off_line = plot_kinetics_disso.line('pl_disso_time_list', 'pl_disso_list', source=kin_source, color='#fc8d62',
                                        name='on_line', line_width=2, line_alpha=0.5)
    circle_kinetics_disso = plot_kinetics_disso.scatter(x='text_x', y='text_y', size=10, color='#fc8d62', source=anno_text_source)
    label_kinetics_disso = LabelSet(x='text_x', y='text_y', x_units='screen', text='disso_text', source=anno_text_source,
                         border_line_color=None, x_offset=0, y_offset=5, border_line_alpha=1.0,
                         background_fill_color=None, background_fill_alpha=1.0)
    plot_kinetics_disso.add_layout(label_kinetics_disso)
    # Widget and callback
    kin_div_kd = Div(
        text="<p>Equilibrium Dissociation Constant (<b><i>K</i><sub>d</sub></b> = <b><i>k</i><sub>off</sub></b>"
             + "/<b><i>k</i><sub>on</sub></b>):<br>" + "{:.5f}".format(kd)
             + " nM</p>", width=300, margin=(5, 0, 5, 5))
    kin_cb = CustomJS(
        args=dict(kin_source=kin_source, hline_source=hline_source, anno_text_source=anno_text_source,
                  shown_points_num=shown_points_num,
                  num_array=num_array, disso_num_array=disso_num_array,
                  kin_div_kd=kin_div_kd, yr_off=plot_kinetics_disso.y_range, yr=plot_kinetics_on.y_range,
                  kin_data_div=kin_data_div, div_text_source=div_text_source),
        code="""
        const kin_source_data = kin_source.data;
        const hline_source_data = hline_source.data;
        const anno_text_source_data = anno_text_source.data;
        const div_text_source_data = div_text_source.data;
        const p0 = kin_source_data['p0'];
        const l0 = kin_source_data['l0'];
        const koff = kin_source_data['koff'];
        const kon = kin_source_data['kon'];
        const time_list = kin_source_data['time_list'];
        const pl_asso_list = kin_source_data['pl_asso_list'];
        const pl_asso_frac_list_vs_pl_1 = kin_source_data['pl_asso_frac_list_vs_pl_1'];
        const pl_asso_frac_list_vs_p0 = kin_source_data['pl_asso_frac_list_vs_p0'];
        const pseudo_pl_asso_list = kin_source_data['pseudo_pl_asso_list'];
        const hline_x = hline_source_data['hline_x'];
        const hline_pl_1 = hline_source_data['hline_pl_1'];
        const hline_pl_eq_pseudo = hline_source_data['hline_pl_eq_pseudo'];
        const text_y = anno_text_source_data['text_y'];
        const text_y_pseudo = anno_text_source_data['text_y_pseudo'];
        const asso_text = anno_text_source_data['asso_text'];
        const asso_text_pseudo = anno_text_source_data['asso_text_pseudo'];
        const disso_text = anno_text_source_data['disso_text'];
        // Dissociation parameters
        const pl_disso_list = kin_source_data['pl_disso_list'];
        const pl_disso_time_list = kin_source_data['pl_disso_time_list'];
        const user_value = cb_obj.value;
        const on_value = user_value * 1e-9;
        if (cb_obj.name == 'kin_spinner_kon') {
            for (var i=0; i<shown_points_num; i++) {
                kon[i] = on_value;
            }
            const kd_html = koff[0] / kon[0];
            kin_div_kd.text="<p>Equilibrium Dissociation Constant (<b><i>K</i><sub>d</sub></b> = "
                            + "<b><i>k</i><sub>off</sub></b>" + "/<b><i>k</i><sub>on</sub></b>):<br>" 
                            + kd_html.toFixed(5) + " nM</p>";
        }
        if (cb_obj.name == 'kin_spinner_koff') {
            for (var i=0; i<shown_points_num; i++) {
                koff[i] = user_value;
            }
            const kd_html = koff[0] / kon[0];
            kin_div_kd.text="<p>Equilibrium Dissociation Constant (<b><i>K</i><sub>d</sub></b> = " 
                            + "<b><i>k</i><sub>off</sub></b>" + "/<b><i>k</i><sub>on</sub></b>):<br>" 
                            + kd_html.toFixed(5) + " nM</p>";
        }
        if (cb_obj.name == 'kin_spinner_p0') {
            for (var i=0; i<shown_points_num; i++) {
                p0[i] = user_value;
            }
        }
        if (cb_obj.name == 'kin_spinner_l0') {
            for (var i=0; i<shown_points_num; i++) {
                l0[i] = user_value;
            }
        }
        // Range update
        yr.start = -0.05 * p0[0];
        yr.end = 1.25 * p0[0];
        const a = 1;
        const b = 0 - (p0[0] + l0[0] + koff[0] / kon[0]);
        const c = p0[0] * l0[0];
        const pl_1 = (0 - b - Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        const pl_2 = (0 - b + Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        const A = 1 / (pl_1 - pl_2);
        const pl_99perc = 0.99 * pl_1;
        const t99 = (A / kon[0]) * Math.log((pl_2 * (pl_1 - pl_99perc)) / (pl_1 * (pl_2 - pl_99perc)));
        for (var i=0; i< shown_points_num; i++){
            pl_asso_list[i] = pl_1 * num_array[i];
            time_list[i] = (A / kon[0]) * Math.log((pl_2 * (pl_1 - pl_asso_list[i])) / (pl_1 * (pl_2 - pl_asso_list[i])));
            pl_asso_frac_list_vs_pl_1[i] = pl_asso_list[i] / pl_1;
            pl_asso_frac_list_vs_p0[i] = pl_asso_list[i] / p0[0];
        }
        // Pseudo rate section
        const pl_eq_pseudo = kon[0] * l0[0] * p0[0] / (kon[0] * l0[0] + koff[0]);
        const kobs = kon[0] * l0[0] + koff[0];
        const pseudo_t99 = Math.log(100)/kobs;
        for (var i=0; i<shown_points_num; i++){
            pseudo_pl_asso_list[i] = (kon[0] * p0[0] * l0[0] / kobs) * (1 - Math.exp(0 - kobs * time_list[i]));
        }
        // Dissociation section
        yr_off.start = -0.05 * pl_1;
        yr_off.end = 1.1 * pl_1;
        const disso_t99 = Math.log(100) / koff[0];
        const disso_t50 = Math.log(2) / koff[0];
        for (var i=0; i<shown_points_num; i++){
            pl_disso_list[i] = pl_1 * disso_num_array[i];
            pl_disso_time_list[i] = (-1 / koff[0] * Math.log(pl_disso_list[i] / pl_1));
        }
        hline_x[1] = t99 * 1.2;
        hline_pl_1[0] = pl_1;
        hline_pl_1[1] = pl_1;
        hline_pl_eq_pseudo[0] = pl_eq_pseudo;
        hline_pl_eq_pseudo[1] = pl_eq_pseudo;
        text_y[0] = pl_1;
        text_y_pseudo[0] = pl_eq_pseudo;
        asso_text[0] = '[PL]eq in 2nd-order process: ' + pl_1.toExponential(2) + ' nM';
        asso_text_pseudo[0] = '[PL]eq in pseudo-1st-order process: ' + pl_eq_pseudo.toExponential(2) + ' nM';
        disso_text[0] = '[PL]0 = ' + pl_1.toExponential(2) + ' nM';
        // Text update
        div_text_source_data['t99'][0] = t99;
        div_text_source_data['pseudo_t99'][0] = pseudo_t99;
        div_text_source_data['pt99_ratio'][0] = pseudo_t99/t99;
        div_text_source_data['disso_t50'][0] = disso_t50;
        div_text_source_data['disso_t99'][0] = disso_t99;
        kin_data_div.text='Time to Reach 99% Equilibrium in Pseudo-1st-order Process ' 
            + '(<b><i>t</i><sub>pseudo-0.99</sub></b>) = ' + div_text_source_data['pseudo_t99'][0].toExponential(2) 
            + ' s;<br>' + 'Time to Reach 99% Equilibrium in 2nd-order Process ' + '(<b><i>t</i><sub>0.99</sub></b>) = ' 
            + div_text_source_data['t99'][0].toExponential(2) + ' s;<br>' 
            + '<b><i>t</i><sub>pseudo-0.99</sub></b>/<b><i>t</i><sub>0.99</sub></b>' 
            + ' = ' + div_text_source_data['pt99_ratio'][0].toFixed(2) + ';<br>' 
            + 'Half-life of Dissociation (ln(2)/koff) = ' 
            + div_text_source_data['disso_t50'][0].toExponential(2) + ' s;<br>' 
            + 'Time to Reach 99% Dissociation = ' + div_text_source_data['disso_t99'][0].toExponential(2) + ' s.';
        kin_source.change.emit();
        hline_source.change.emit();
        anno_text_source.change.emit();
        div_text_source.change.emit();
    """)
    # Configuration parameters and widgets
    # In the webpage, the unit of kon is M-1 s-1; in backend, the unit of kon is nM-1 s-1.
    kin_header = Div(text='<h1 style="display:inline; margin-top: 0">Kinetics of Association and '
                          +'Dissociation</h1>'
                          + '<span>&nbsp;<a href="/binding-curve-viewer/">[Back to Home]</a></span>',
                     height=40, width=800, height_policy='auto', margin=(5, 5, 5, 5), )
    kin_div_p0 = Div(text="<p>Initial Concentration of Fixed Component (e.g., Protein, <b>[P]<sub>0</sub></b>, nM): "
                         + "<br>0.001 ~ 50000 nM</p>", width=300, height=48, margin=(0, 0, 0, 5))
    kin_spinner_p0 = Spinner(low=1e-3, high=5e4, value=p0, step=1, title='', name='kin_spinner_p0', width=240)
    kin_spinner_p0.js_on_change('value', kin_cb)
    kin_div_l0 = Div(text="<p>Initial Concentration of Ligand (<b>[L]<sub>0</sub></b>, nM):" +
                          "<br>0.001 ~ 50000 nM</p>", width=300, height=30, margin=(5, 0, 0, 5))
    kin_spinner_l0 = Spinner(low=1e-3, high=5e5, value=l0, step=1, title='', name='kin_spinner_l0', width=240)
    kin_spinner_l0.js_on_change('value', kin_cb)
    kin_div_koff = Div(text="<p>Off-rate Constant (<b><i>k</i><sub>off</sub></b>, s<sup>-1</sup>): "
                           + "1e-6 ~ 1 s<sup>-1</sup></p>", width=300, height=15, margin=(5, 0, 0, 5))
    kin_spinner_koff = Spinner(low=1e-6, high=1, value=koff, step=1e-6, title='',
                              name='kin_spinner_koff', width=240)
    kin_spinner_koff.js_on_change('value', kin_cb)
    kin_div_kon = Div(text="<p>On-rate Constant (<b><i>k</i><sub>on</sub></b>, M<sup>-1</sup> s<sup>-1</sup>): "
                           + "1e3 ~ 1e8 M<sup>-1</sup> s<sup>-1</sup></p>", width=300, height=15,
                      margin=(5, 0, 0, 5))
    kin_spinner_kon = Spinner(low=1e3, high=1e8, value=kon * 1e9, step=1e3, title='',
                              name='kin_spinner_kon', width=240)
    kin_spinner_kon.js_on_change('value', kin_cb)
    kin_config_col = column(kin_div_p0, kin_spinner_p0, kin_div_l0, kin_spinner_l0,
                           kin_div_koff, kin_spinner_koff, kin_div_kon, kin_spinner_kon,
                           kin_div_kd, kin_data_div)
    return (kin_header, kin_config_col, plot_kinetics_on, plot_kinetics_disso)


def iso_affinity_plot(p0=100, l0=150, plot_points=100, shown_points_num=100):
    r"""To further compare different kinetic properties under the condition of :math:`[\mathrm{P}]_0`
    = 100 nM and :math:`[\mathrm{L}]_0` = 150 nM, we used this function to generate the iso-affinity
    graph of the second-order binding process. The resultant HTML file will show three parts:
    1. iso-affinity graph; 2. binding curves; 3. a table of experiment information.

    Parameters
    ----------
    p0 : float
        The initial concentration of the free protein in the unit of nM (default value is 100).
    l0 : float
        The initial concentration of the free ligand in the unit of nM (default value is 150).
    plot_points : int
        The number of points (i.e., smoothness) of lines (default value is 100).
    shown_points_num : int
        The number of points in lines shown in the figure (default value is 100). To avoid the
        long tail before equilibrium, the binding curve ends when 99.0% equilibrium is reached,
        shown_points_num should be less than plot_points.

    Returns
    -------
    Tuple
        A tuple of widgets and figures, including ``p3_header``, ``fig_points``,
        ``fig_curves``, ``kin_div``, ``kin_table`` can be used to generate the web application.
    """
    # The order of data point preparation (x, y)
    # (1,1 (2,1) (3,1) ... (6,1)
    # -> (1,2) (2,2) (3,2) ... (6,2)
    # -> (1,3) ... (6,3)
    # -> ... -> (1,7) ... (6,7)
    kon_list = []
    kon_list.extend([1e3, 1e4, 1e5, 1e6, 1e7, 1e8] * 7)
    koff_list = []
    koff_list.extend([1e-6] * 6)
    koff_list.extend([1e-5] * 6)
    koff_list.extend([1e-4] * 6)
    koff_list.extend([1e-3] * 6)
    koff_list.extend([1e-2] * 6)
    koff_list.extend([1e-1] * 6)
    koff_list.extend([1] * 6)
    kd_list = [
        1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14,
        1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13,
        1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12,
        1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11,
        1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10,
        1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9,
        1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
    data_points = {'x': kon_list, 'y': koff_list, 'z': kd_list}
    source_points = ColumnDataSource(data_points)
    # Point figure generation
    fig_points = figure(width=550, height=450, x_axis_type='log', y_axis_type='log', tools=[], output_backend='svg')
    fig_points.toolbar.logo = None
    fig_points.grid.visible = False
    fig_points.outline_line_alpha = 0
    fig_points.xaxis.axis_label = r"$$k_{\text{on}} (\text{M}^{-1} \text{s}^{-1})$$"
    fig_points.yaxis.axis_label = r"$$k_{\text{off}} (\text{s}^{-1})$$"
    fig_points.xaxis.bounds = (1e3, 1e8)
    fig_points.yaxis.bounds = (1e-6, 1)
    # Main points
    points_res = fig_points.scatter(x='x', y='y', source=source_points, size=25, nonselection_alpha=0.5,
                                   fill_color='gray', fill_alpha=0.5, line_width=0.1)
    points_res.hover_glyph = Scatter(fill_color='gray', size=25, line_width=2, line_color='orange')
    # Helper circle for display range definition
    fig_points.scatter(x=(7e8), y=(8e-7), size=1, fill_alpha=0, line_color='white')
    # Layer generation
    color_list = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f'] * 6
    layer_renderer = []
    for i in range(len(kon_list)):
        layer_renderer.append(
            fig_points.scatter(x=kon_list[i], y=koff_list[i], size=25, fill_color=color_list[i], line_width=0.1)
        )
    for p in layer_renderer:
        p.visible = True
    # Curve figure generation
    line_width_list = []
    line_width_list.extend([6.8] * 6)
    line_width_list.extend([6.0] * 6)
    line_width_list.extend([5.2] * 6)
    line_width_list.extend([4.4] * 6)
    line_width_list.extend([3.6] * 6)
    line_width_list.extend([2.8] * 6)
    line_width_list.extend([2.0] * 6)
    table_data = dict(koff=[], kon=[], kd=[], pl_1=[], t50=[], t99=[], point_color=[])
    fig_curves = figure(width=450, height=450, x_axis_type='log', tools=[], y_range=(-0.04, 1), output_backend='svg')
    fig_curves.toolbar.logo = None
    fig_curves.xaxis.axis_label = 'Time (s)'
    fig_curves.yaxis.axis_label = r"Fraction of Bound Protein ($$\text{[PL]}/\text{[P]}_0$$)"
    fig_curves.xaxis.axis_label_text_font_style = 'normal'
    fig_curves.yaxis.axis_label_text_font_style = 'normal'
    curve_renderer = []
    for i in range(len(kon_list)):
        kon = data_points['x'][i] * 1e-9  # nM-1 s-1
        koff = data_points['y'][i]
        kd = data_points['z'][i] * 1e9  # nM
        # a, b, c are the coefficients of quadratic equation
        # For more details, see the supplementary information of binding curve viewer paper
        a = 1
        b = 0 - (p0 + l0 + koff / kon)
        c = p0 * l0
        pl_1 = (0 - b - np.sqrt(b * b - 4 * a * c)) / (2 * a)  # pl_1 is [PL]eq, pl_1 < pl_2
        pl_2 = (0 - b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
        A = 1 / (pl_1 - pl_2)
        pl_50perc = pl_1 * 0.50  # [PL] @ half equilibrium
        pl_99perc = pl_1 * 0.99  # [PL] @ 99% equilibrium
        thalf =  A / kon * np.log((pl_2 * (pl_1 - pl_50perc)) / (pl_1 * (pl_2 - pl_50perc)))  # time for 50% equilibrium
        t99 = A / kon * np.log((pl_2 * (pl_1 - pl_99perc)) / (pl_1 * (pl_2 - pl_99perc)))  # time for 99% equilibrium
        table_data['koff'].append('%.2e' % data_points['y'][i])
        table_data['kon'].append('%.2e' % data_points['x'][i])
        table_data['kd'].append('%.2e' % data_points['z'][i])
        table_data['pl_1'].append('%.2e' % pl_1)
        table_data['t99'].append('%.2e' % t99)
        table_data['t50'].append('%.2e' % thalf)
        table_data['point_color'].append(color_list[i])
        num_array = np.arange(0, plot_points, 1) / plot_points
        pl_asso_list = pl_1 * num_array[:shown_points_num]
        time_list = A / kon * np.log((pl_2 * (pl_1 - pl_asso_list)) / (pl_1 * (pl_2 - pl_asso_list)))
        frac_pl_list = pl_asso_list / p0

        tmp_data = dict(t=time_list, y=frac_pl_list, kon=[data_points['x'][i]] * shown_points_num,
                        koff=[data_points['y'][i]] * shown_points_num, kd=[data_points['z'][i]] * shown_points_num,
                        t99=[t99] * shown_points_num, f99=[frac_pl_list[-1]] * shown_points_num, in_index=[i] * shown_points_num)
        tmp_source = ColumnDataSource(tmp_data)
        curve_renderer.append(
            fig_curves.line(x='t', y='y', source=tmp_source, line_width=line_width_list[i], line_color=color_list[i])
        )
    for p in curve_renderer:
        p.visible = True
    fig_curves.add_tools(HoverTool(
        tooltips=[('Index', '@in_index'),
                  ('koff (s-1)', '@koff{%0.1e}'),
                  ('kon (M-1 s-1)', '@kon{%0.1e}'),
                  ('Kd (M)', '@kd{%0.1e}'),
                  ('Time to 99% equilibrium (s)', '@t99{%0.1e}'),
                  ('[P]/[P]0 at 99% equilibrium', '@f99{0.2e}')],
        formatters={'@koff': 'printf', '@kon': 'printf', '@kd': 'printf', '@t99': 'printf'},
    ))
    # Annotation dash line and label
    fig_points.line((5e2, 2e4), (5e-2, 2), line_dash='dotdash', line_color='gray', line_width=3)
    fig_points.line((5e2, 2e5), (5e-3, 2), line_dash='dotdash', line_color='gray', line_width=3)
    fig_points.line((5e2, 2e6), (5e-4, 2), line_dash='dotdash', line_color='gray', line_width=3)
    fig_points.line((5e2, 2e7), (5e-5, 2), line_dash='dotdash', line_color='gray', line_width=3)
    fig_points.line((5e2, 2e8), (5e-6, 2), line_dash='dotdash', line_color='gray', line_width=3)
    fig_points.line((5e2, 2e8), (5e-7, 2e-1), line_dash='dotdash', line_color='gray', line_width=3)
    fig_points.line((5e3, 2e8), (5e-7, 2e-2), line_dash='dotdash', line_color='gray', line_width=3)
    fig_points.line((5e4, 2e8), (5e-7, 2e-3), line_dash='dotdash', line_color='gray', line_width=3)
    fig_points.line((5e5, 2e8), (5e-7, 2e-4), line_dash='dotdash', line_color='gray', line_width=3)
    fig_points.line((5e6, 2e8), (5e-7, 2e-5), line_dash='dotdash', line_color='gray', line_width=3)
    label_source = ColumnDataSource(
        data=dict(x=[600, 2e3, 2e4, 2e5, 2e6, 2e7, 2e8, 2e8, 2e8, 2e8, 2e8, 2e8, 2e8],
                  y=[2, 2, 2, 2, 2, 2, 2, 2e-1, 2e-2, 2e-3, 2e-4, 2e-5, 1e-6],
                  text=['Kd', '1 mM', '100 uM', '10 uM', '1 uM', '100 nM',
                        '10 nM', '1 nM', '100 pM', '10 pM', '1 pM', '100 fM', '10 fM']))
    labels = LabelSet(x='x', y='y', text='text', text_color='gray',
                      x_offset=-9, y_offset=2, source=label_source)
    fig_points.add_layout(labels)
    # Table widget and callback
    table_source = ColumnDataSource(table_data)
    table_view = CDSView(filter=BooleanFilter(booleans=[True] * 42))
    table_columns = [
        TableColumn(field='kd', title='Kd (M)',
                    formatter=HTMLTemplateFormatter(
                        template='<div style="background:<%= point_color %>"><%= value %></div>')),
        TableColumn(field='koff', title='koff (s-1)'),
        TableColumn(field='kon', title='kon (M-1 s-1)'),
        TableColumn(field='pl_1', title='[PL]eq (nM)'),
        TableColumn(field='t50', title='t0.50 (s)'),
        TableColumn(field='t99', title='t0.99 (s)'), ]
    kin_table = DataTable(columns=table_columns, source=table_source, view=table_view,
                          width=480, height=350, reorderable=False, sortable=False)
    # Divider and Kinetics description
    p3_header = Div(text='<h1 style="display:inline; margin-top: 0">Iso-affinity Graph</h1>'
                    + '<span>&nbsp;<a href="/binding-curve-viewer/">[Back to Home]</a></span>',
                    height=45)
    kin_div = Div(text='The time to 50% and 99% equilibrium and <b>[PL]<sub>eq</sub></b> '
                       + 'under specific initial conditions (<b>[P]<sub>0</sub></b> = 100 nM '
                       + 'and <b>[L]<sub>0</sub></b> = 150 nM) and different kinetic parameters '
                       + 'in the second-order binding process:',
                  width=480, height=55)
    # Code and event
    code_sele = '''
    if (source_points.selected.indices.length > 0){
        var idx = source_points.selected.indices;
        var selected_index = idx[0];
        var b = lines[selected_index].visible;
        var fb = view.filter.booleans[selected_index];
        lines[selected_index].visible = ! b;
        layer[selected_index].visible = ! b;
        view.filter.booleans[selected_index] = ! fb;
        view.filter.change.emit();
    }'''
    code_clear = '''
    if (source_points.selected.indices.length == 0){
        for (var i=0; i < lines.length; i++){
            lines[i].visible = false;
            layer[i].visible = false;
            view.filter.booleans[i] = false;
        }
        view.filter.change.emit();
    }'''

    fig_points.js_on_event('tap', CustomJS(args={'source_points': source_points, 'lines': curve_renderer,
                                                 'view': kin_table.view, 'layer': layer_renderer}, code=code_clear))
    tap_points = TapTool(renderers=[points_res], mode='replace',
                  callback=CustomJS(args={'source_points': source_points, 'lines': curve_renderer,
                                          'view': kin_table.view, 'layer': layer_renderer}, code=code_sele))
    hover = HoverTool(renderers=[points_res],
                      tooltips=[('Index', "$index"),
                                ('koff (s-1)', '@y{%0.1e}'),
                                ('kon (M-1 s-1)', '@x{%0.1e}'),
                                ('Kd (M)', '@z{%0.1e}')],
                      formatters={"@x": "printf", "@y": "printf", "@z": "printf"})
    fig_points.add_tools(tap_points)
    fig_points.add_tools(hover)
    return p3_header, fig_points, fig_curves, kin_div, kin_table


def competitive_binding_association(p0=100, l0=75, i0=75, pl_koff=1e-2, pl_kon=1e5, pi_koff=1e-2, pi_kon=1e5,
                                    diff_calc_resolution=200000, t99Ntimes=10):
    r"""The simulation of the competitive binding of association is initiated by mixing an arbitrary
    volume of protein with an arbitrary volume of the mixture of ligand and inhibitor. The working
    concentration of each component at the time of mixing can be determined by the p0, l0, and i0.

    .. math:: :label: eq_4.1

        \mathrm{P} + \mathrm{L} \mathrel{\mathop{\rightleftarrows}^{k_{\mathrm{ligand\textnormal{-}off}}}_{k_{\mathrm{ligand\textnormal{-}on}}}} \mathrm{PL}

    .. math:: :label: eq_4.2

        \mathrm{P} + \mathrm{I} \mathrel{\mathop{\rightleftarrows}^{k_{\mathrm{inhibitor\textnormal{-}off}}}_{k_{\mathrm{inhibitor\textnormal{-}on}}}} \mathrm{PI}

    The competitive binding is described by above equations. To describe the kinetic process of
    the competitive binding, we used the Mahan-Motulsky equation to calculate the relationship
    of the :math:`[\mathrm{PL}]` and time in the association process. The Mahan-Motulsky equation can only be
    applied in the pseudo-first-order binding process. We devised the numerical method to
    simulate the second-order kinetic process of the competitive binding of both association
    and dissociation. See the Methods of :ref:`the paper<bcv_paper>` for detailed expression.

    Parameters
    ----------
    p0 : float
        The initial concentration of the free protein in the unit of nM (default value is 100).
    l0 : float
        The initial concentration of the free ligand in the unit of nM (default value is 75).
    i0 : float
        The initial concentration of the free inhibitor in the unit of nM (default value is 75).
    pl_koff : float
        The value of the ligand's dissociation rate constant (:math:`k_{\mathrm{ligand\textnormal{-}off}}`) in the unit of s\ :sup:`-1`
        (default value is 1e-2).
    pl_kon : float
        The value of the ligand's association rate constant (:math:`k_{\mathrm{ligand\textnormal{-}on}}`) in the unit of M\ :sup:`-1` s\ :sup:`-1`
        (default value is 1e5).
    pi_koff : float
        The value of the inhibitor's dissociation rate constant (:math:`k_{\mathrm{inhibitor\textnormal{-}off}}`) in the unit of s\ :sup:`-1`
        (default value is 1e-2).
    pi_kon : float
        The value of the inhibitor's association rate constant (:math:`k_{\mathrm{inhibitor\textnormal{-}on}}`) in the unit of M\ :sup:`-1` s\ :sup:`-1`
        (default value is 1e5).
    diff_calc_resolution : int
        The resolution of the numerical simulation (default value is 200,000).
    t99Ntimes : int
        The value affects the time length of the simulation (default value is 10).

    Returns
    -------
    Tuple
        A tuple of widgets and figures, including ``p4_header``, ``pro_div``, ``lig_div``, ``inhib_div``,
        ``pro_spinner``, ``lig_spinner``, ``inhib_spinner``, ``lig_koff_div``, ``lig_koff_spinner``,
        ``lig_kon_div``,  ``lig_kon_spinner``, ``inhib_koff_div``, ``inhib_koff_spinner``,  ``inhib_kon_div``,
        ``inhib_kon_spinner``, ``info_div``, ``error_div``, ``ic50_div``, ``ref_div``,  ``plot_kinetics_comp_asso``,
        ``ic50_fig`` can be used to generate the web application.
    """

    def second_order_kinetics(p0, l0, koff, kon):
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
        return time_list, asso_list, t99, pl_1

    def second_order_kinetics_rerun(p0, l0, koff, kon, t99):
        # This function use the t99 to generate PL list
        kon = kon * 1e-9  # convert the unit from M^-1 s^-1 to nM^-1 s^-1
        a = 1
        b = 0 - (p0 + l0 + koff / kon)
        c = p0 * l0
        pl_1 = (0 - b - np.sqrt(b * b - 4 * a * c)) / (2 * a)  # pl_1 is [PL]eq, pl_1 < pl_2
        pl_2 = (0 - b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
        num_array = np.linspace(0, 1, 100)  # take in the t99 into range
        time_list = t99 * num_array
        C = (pl_1 / pl_2) * np.exp(kon * (pl_1 - pl_2) * time_list)
        asso_list = (C * pl_2 - pl_1) / (C - 1)
        return time_list, asso_list

    time_list_pl, asso_list_pl, t99_pl, pl_eq_second_order = second_order_kinetics(p0, l0, pl_koff, pl_kon)
    time_list_pi, asso_list_pi, t99_pi, pi_eq_second_order = second_order_kinetics(p0, i0, pi_koff, pi_kon)
    t99_short = None
    t99_long = None
    if t99_pl >= t99_pi:
        t99_long = t99_pl
        t99_short = t99_pi
        time_list_pi, asso_list_pi = second_order_kinetics_rerun(p0, i0, pi_koff, pi_kon, t99_pl)
    else:
        t99_long = t99_pi
        t99_short = t99_pl
        time_list_pl, asso_list_pl = second_order_kinetics_rerun(p0, l0, pl_koff, pl_kon, t99_pi)
    second_order_pl = dict(time_list_pl=time_list_pl, asso_list_pl=asso_list_pl)
    second_order_pi = dict(time_list_pi=time_list_pi, asso_list_pi=asso_list_pi)

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

    def diff_calc_association(p0, l0, i0, pl_koff, pl_kon, pi_koff, pi_kon, pl_comp_eq, t99_short, t99_long,
                              resolution=diff_calc_resolution, t99Ntimes=t99Ntimes, plot_points=300):
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
        time_loop = np.append(time_loop1, time_loop2[1:])  # length should be 2 * resolution - 1
        pl_99 = 0.99 * pl_comp_eq
        pl_101 = 1.01 * pl_comp_eq
        idx_99 = None  # the index of approximate pl at 99% equilibrium
        # the length of consumlist_l is 2 * resolution - 2
        for i in range(0, 2 * resolution - 2):
            pl_list.append(consumlist_l[i] - coffsumlist_l[i])
        for i in reversed(range(2 * resolution - 1)):
            if pl_list[i] < pl_99 or pl_list[i] > pl_101:
                idx_99 = i
                break
        t99_approx_idx = min(range(len(time_loop)), key=lambda i: abs(time_loop[i] - t99_long))
        samples_idx = np.linspace(0, t99_approx_idx, plot_points).astype(int)
        samples_time = [time_loop[idx] for idx in samples_idx]
        samples_pl = [pl_list[idx] for idx in samples_idx]
        return samples_time, samples_pl, time_loop[idx_99]


    def kinetics_pseudo_asso(N, L0, I0, k1, k2, k3, k4, t99, plot_points=300):
        k1 = k1 * 1e-9  # kon_L, nM^-1 s^-1
        k3 = k3 * 1e-9  # kon_I, nM^-1 s^-1
        K_A = k1 * L0 + k2
        K_B = k3 * I0 + k4
        tmp_sqrt = np.sqrt((K_A - K_B) ** 2 + 4 * k1 * k3 * L0 * I0)
        K_F = 0.5 * (K_A + K_B + tmp_sqrt)
        K_S = 0.5 * (K_A + K_B - tmp_sqrt)
        PL_eq = N * k1 * k4 * L0 / (K_F * K_S)
        PL_99 = 0.99 * PL_eq
        PL_101 = 1.01 * PL_eq
        t_N = 2
        t99_2 = t99 * t_N
        comp_tlist = np.linspace(0, t99_2, plot_points * t_N)
        PL = (N * k1 * L0 / (K_F - K_S)) * (
                k4 * (K_F - K_S) / (K_F * K_S) + ((k4 - K_F) / K_F) * np.exp(0 - K_F * comp_tlist)
                - ((k4 - K_S) / K_S) * np.exp(0 - K_S * comp_tlist))
        idx_99 = None
        for i in reversed(range(len(PL))):
            if PL[i] < PL_99 or PL[i] > PL_101:
                idx_99 = i
                break
        return comp_tlist[:plot_points], PL[:plot_points], comp_tlist[idx_99], PL_eq

    def ic50_curve(p0, l0, pl_kd, pi_kd):
        # The units of p0 and l0 are nM
        # The units of pl_kd and pi_kd are nM
        pl_kd_M, pi_kd_M = pl_kd * 1e-9, pi_kd * 1e-9
        # without inhibitor the PL_eq equals pl_1
        a = 1
        b = 0 - (p0 + l0 + pl_kd)
        c = p0 * l0
        pl_1 = (0 - b - np.sqrt(b * b - 4 * a * c)) / (2 * a)
        ic50_approx = pi_kd * (1 + l0 / pl_kd)
        if ic50_approx > 1e6:
            return None, None, None, None, None
        plot_ic50_i_log_list = np.linspace(-4, np.log10(ic50_approx * 100), 100)
        plot_ic50_i_list = []
        plot_ic50_i_free_list = []
        plot_ic50_pl_list = []
        for i_log in plot_ic50_i_log_list:
            plot_ic50_i_list.append(10 ** i_log)
        for i in plot_ic50_i_list:
            p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(p0, l0, i, pl_kd_M, pi_kd_M)
            plot_ic50_pl_list.append(pl_comp_eq)
            plot_ic50_i_free_list.append(i - pi_eq)
        ic50_approx = pi_kd * (1 + l0 / pl_kd)
        inner_ic50_i_list = np.linspace(0, 100 * ic50_approx, 90000)  # this is I total
        inner_ic50_i_free_list = []  #this is I free
        inner_ic50_pl_list = []
        for i in inner_ic50_i_list:
            p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(p0, l0, i, pl_kd_M, pi_kd_M)
            # pl_quan = p0 * l0 / (kd * (1 + i / ki) + l0)
            inner_ic50_pl_list.append(pl_comp_eq)
            inner_ic50_i_free_list.append(i - pi_eq)
        ic50_approx_idx = min(range(len(inner_ic50_pl_list)), key=lambda i: abs(inner_ic50_pl_list[i] - (pl_1 / 2)))
        apparent_ic50 = inner_ic50_i_list[ic50_approx_idx]
        ic50 = inner_ic50_i_free_list[ic50_approx_idx]
        return plot_ic50_i_list, plot_ic50_i_free_list, plot_ic50_pl_list, apparent_ic50, ic50

    p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(P0=p0, A0=l0, B0=i0, KA=pl_koff / pl_kon, KB=pi_koff / pi_kon)
    samples_time, samples_pl, t99_diff = diff_calc_association(p0=p0, l0=l0, i0=i0, pl_koff=pl_koff, pl_kon=pl_kon,
                                                               pi_koff=pi_koff, pi_kon=pi_kon, pl_comp_eq=pl_comp_eq,
                                                               t99_short=t99_short, t99_long=t99_long)
    comp_tlist, PL, t99_pseudo, pl_comp_eq_pseudo = kinetics_pseudo_asso(N=p0, L0=l0, I0=i0, k1=pl_kon, k2=pl_koff,
                                                                         k3=pi_kon, k4=pi_koff, t99=t99_long)
    diff_pl = dict(samples_time=samples_time, samples_pl=samples_pl)
    pseudo_pl = dict(comp_tlist=comp_tlist, PL=PL)
    pl_kd = pl_koff / pl_kon * 1e9
    pi_kd = pi_koff / pi_kon * 1e9
    plot_ic50_i_list, plot_ic50_i_free_list, plot_ic50_pl_list, apparent_ic50, ic50 = ic50_curve(p0, l0, pl_kd, pi_kd)

    def calc_ki_app4(pt, lt, it, pl_kd):
        a = 1
        b = 0 - (pt + lt + pl_kd)
        c = pt * lt
        pl_1 = (0 - b - np.sqrt(b * b - 4 * a * c)) / (2 * a)
        # use the solutions in Anal. Biochem. 2004, 332 (2), 261-73.
        # pt - pl_1 is the P0 in original paper.
        # pl_1 is the PL0 in original paper.
        P0 = pt - pl_1
        PL50 = pl_1 / 2
        L50 = lt - PL50
        I50 = it - pt + pl_kd * PL50 / L50 + PL50
        ki = I50 / (L50 / pl_kd + P0 / pl_kd + 1)
        return ki
    calc_ki = calc_ki_app4(p0, l0, apparent_ic50, pl_kd)

    ic50_data = dict(plot_ic50_i_list=plot_ic50_i_list, plot_ic50_i_free_list=plot_ic50_i_free_list,
                     plot_ic50_pl_list=plot_ic50_pl_list)
    info_data = dict(p0=[p0], l0=[l0], i0=[i0], pl_kd=[pl_kd], pi_kd=[pi_kd], ic50=[ic50], apparent_ic50=[apparent_ic50],
                     pl_koff=[pl_koff], pl_kon=[pl_kon], pi_koff=[pi_koff], pi_kon=[pi_kon],
                     pl_comp_eq=[pl_comp_eq], pl_comp_eq_pseudo=[pl_comp_eq_pseudo],
                     pl_eq_second_order=[pl_eq_second_order], pi_eq_second_order=[pi_eq_second_order],
                     t99_short=[t99_short], t99_long=[t99_long],
                     t99_pl=[t99_pl], t99_pi=[t99_pi], t99_diff=[t99_diff], t99_pseudo=[t99_pseudo],
                     diff_calc_resolution=[diff_calc_resolution], t99Ntimes=[t99Ntimes])
    second_order_pl_source = ColumnDataSource(data=second_order_pl)
    second_order_pi_source = ColumnDataSource(data=second_order_pi)
    diff_pl_source = ColumnDataSource(data=diff_pl)
    pseudo_pl_source = ColumnDataSource(data=pseudo_pl)
    ic50_source = ColumnDataSource(data=ic50_data)
    info_source = ColumnDataSource(data=info_data)
    if pl_eq_second_order > pi_eq_second_order:
        pl_eq_large = pl_eq_second_order
    else:
        pl_eq_large = pi_eq_second_order
    plot_kinetics_comp_asso = figure(width=480, height=575, tools=[], output_backend='svg',
                                     y_range=(-0.05 * pl_eq_large, 1.5 * pl_eq_large), )
    plot_kinetics_comp_asso.toolbar.logo = None
    plot_kinetics_comp_asso.xaxis.axis_label = "Time (s)"
    plot_kinetics_comp_asso.xaxis.axis_label_text_font_style = 'normal'
    plot_kinetics_comp_asso.yaxis.axis_label_text_font_style = 'normal'
    plot_kinetics_comp_asso.yaxis.axis_label = "Concentration of Bound Protein (nM)"
    plot_kinetics_comp_asso.axis.minor_tick_in = -3
    plot_kinetics_comp_asso.axis.minor_tick_out = 6
    pro_lig_line = plot_kinetics_comp_asso.line('time_list_pl', 'asso_list_pl', source=second_order_pl_source,
                                                color='#fc8d62', line_width=2, line_alpha=0.8, line_dash='dashed',
                                                legend_label='Only ligand, [PL]')
    pro_inhib_line = plot_kinetics_comp_asso.line('time_list_pi', 'asso_list_pi', source=second_order_pi_source,
                                                  color='black', line_width=2, line_alpha=0.8, line_dash='dotdash',
                                                  legend_label='Only inhibitor, [PI]')
    comp_line_diff = plot_kinetics_comp_asso.line('samples_time', 'samples_pl', source=diff_pl_source,
                                                  color='#fc8d62', line_width=3, line_alpha=0.8,
                                                  legend_label='Competition-simulation, [PL]')
    comp_line_pseudo = plot_kinetics_comp_asso.line('comp_tlist', 'PL', source=pseudo_pl_source,
                                                    color='#66c2a5', line_width=3, line_alpha=0.8,
                                                    legend_label='Competition-pseudo, [PL]')
    plot_kinetics_comp_asso.legend.location = 'top_left'


    ic50_fig = figure(width=480, height=452, x_axis_type='log', tools=[], output_backend='svg')
    ic50_fig.toolbar.logo = None
    ic50_fig.xaxis.axis_label = "Concentration of Free/Total Inhibitor (nM)"
    ic50_fig.xaxis.axis_label_text_font_style = 'normal'
    ic50_fig.yaxis.axis_label_text_font_style = 'normal'
    ic50_fig.yaxis.axis_label = "Concentration of Protein-Ligand Complex ([PL], nM)"
    ic50_fig.axis.minor_tick_in = -3
    ic50_fig.axis.minor_tick_out = 6
    apparent_ic50_line = ic50_fig.line('plot_ic50_i_list', 'plot_ic50_pl_list', source=ic50_source,
                              color='#fc8d62', line_width=2, line_alpha=0.8, legend_label='[I] total')
    ic50_line = ic50_fig.line('plot_ic50_i_free_list', 'plot_ic50_pl_list', source=ic50_source,
                              color='#66c2a5', line_width=2, line_alpha=0.8, legend_label='[I] free')
    ic50_fig.legend.location = 'bottom_left'
    info_div = Div(
        text=
        '<b><i>K</i><sub>d</sub></b> for Protein and Ligand (<b><i>K</i><sub>d</sub></b> = <b><i>k</i><sub>off-ligand</sub></b>/<b><i>k</i><sub>on-ligand</sub></b>): ' + '{:.5f}'.format(
            pl_kd) + ' nM;<br>' +
        '<b><i>K</i><sub>d</sub></b> for Protein and Inhibitor (<b><i>K</i><sub>d</sub></b>, or <b><i>K</i><sub>i</sub></b> = <b><i>k</i><sub>off-inhibitor</sub></b>/<b><i>k</i><sub>on-inhibitor</sub></b>): ' + '{:.5f}'.format(
            pi_kd) + ' nM;<br>' +
        'Exact Time to Reach 99% Equilibrium (<b><i>t</i><sub>0.99</sub></b>), <b>[PL]<sub>eq</sub></b>, and <b>[PI]<sub>eq</sub></b> in the Second-order Association Process of Only Ligand or Inhibitor:<br>' +
        '<b><i>t</i><sub>0.99-ligand</sub></b> = ' + '%.2e' % t99_pl + ' s, <b><i>t</i><sub>0.99-inhibitor</sub></b> = ' + '%.2e' % t99_pi + ' s,<br>' +
        '<b>[PL]<sub>eq</sub></b> = ' + '{:.5f}'.format(
            pl_eq_second_order) + ' nM, <b>[PI]<sub>eq</sub></b> = ' + '{:.5f}'.format(
            pi_eq_second_order) + ' nM;<br>' +
        'Approximate <b><i>t</i><sub>0.99</sub></b> in Competitive Second-order Association Process: ' + '%.2e' % t99_diff + ' s;<br>' +
        'Approximate <b><i>t</i><sub>0.99</sub></b> in Competitive Pseudo-first-order Association Process: ' + '%.2e' % t99_pseudo + ' s;<br>' +
        'Exact <b>[PL]<sub>eq</sub></b> in Competitive Second-order Binding Process: ' + '{:.5f}'.format(pl_comp_eq) + ' nM;<br>' +
        'Exact <b>[PL]<sub>eq</sub></b> in Competitive Pseudo-first-order Binding Process: ' + '{:.5f}'.format(
            pl_comp_eq_pseudo) + ' nM;<br>',
        width=500, margin=(15, 25, 5, 5))
    error_div = Div(text='', width=500)
    ic50_div = Div(
        text=
        'Under the Condition of <b>[P]<sub>0</sub></b> = ' + '{:.5f}'.format(
            p0) + ' nM, <b>[L]<sub>0</sub></b> = ' + '{:.5f}'.format(l0) +
        ' nM, <b><i>K</i><sub>d</sub></b> = ' + '{:.5f}'.format(pl_kd) +
        ' nM, and <b><i>K</i><sub>i</sub></b> = ' + '{:.5f}'.format(
            pi_kd) + ' nM, <b>IC<sub>50</sub></b> Approximately Equals ' + '{:.5f}'.format(ic50) +
        ' nM, Apparent <b>IC<sub>50</sub></b> Approximately Equals ' + '{:.5f}'.format(apparent_ic50) +
        ' nM by Wang Equation<sup>1</sup>;<br>' +'<b><i>K</i><sub>i</sub></b> Calculated by Equations<sup>2</sup> Equals '+'{:.5f} nM.'.format(calc_ki),
        width=440, margin=(5, 5, 13, 35))
    ref_div = Div(
        text='<b>References</b><br>' + '[1] FEBS Lett. 2003, 360 (2), 111-114.<br>[2] Anal. Biochem. 2004, 332 (2), 261-273.'
    )
    comp_callback = CustomJS(args=dict(second_order_pl_source=second_order_pl_source,
                                       second_order_pi_source=second_order_pi_source,
                                       diff_pl_source=diff_pl_source, pseudo_pl_source=pseudo_pl_source,
                                       ic50_source=ic50_source,
                                       info_source=info_source, info_div=info_div,
                                       error_div=error_div, ic50_div=ic50_div,
                                       yr=plot_kinetics_comp_asso.y_range, ), code="""
    const second_order_pl_source_data = second_order_pl_source.data;
    const time_list_pl = second_order_pl_source_data['time_list_pl'];
    const asso_list_pl = second_order_pl_source_data['asso_list_pl'];
    const second_order_pi_source_data = second_order_pi_source.data;
    const time_list_pi = second_order_pi_source_data['time_list_pi'];
    const asso_list_pi = second_order_pi_source_data['asso_list_pi'];
    const diff_pl_source_data = diff_pl_source.data;
    const samples_time = diff_pl_source_data['samples_time'];
    const samples_pl = diff_pl_source_data['samples_pl'];
    const pseudo_pl_source_data = pseudo_pl_source.data;
    const comp_tlist = pseudo_pl_source_data['comp_tlist'];
    const PL = pseudo_pl_source_data['PL'];
    const ic50_source_data = ic50_source.data;
    const plot_ic50_i_list = ic50_source_data['plot_ic50_i_list'];
    const plot_ic50_i_free_list = ic50_source_data['plot_ic50_i_free_list'];
    const plot_ic50_pl_list = ic50_source_data['plot_ic50_pl_list'];

    const info_source_data = info_source.data;
    const p0 = info_source_data['p0'];
    const l0 = info_source_data['l0'];
    const i0 = info_source_data['i0'];
    const pl_kd = info_source_data['pl_kd'];
    const pi_kd = info_source_data['pi_kd'];
    const ic50 = info_source_data['ic50'];
    const pl_koff = info_source_data['pl_koff'];
    const pl_kon = info_source_data['pl_kon'];
    const pi_koff = info_source_data['pi_koff'];
    const pi_kon = info_source_data['pi_kon'];
    const pl_comp_eq = info_source_data['pl_comp_eq'];
    const pl_comp_eq_pseudo = info_source_data['pl_comp_eq_pseudo'];
    const pl_eq_second_order = info_source_data['pl_eq_second_order'];
    const pi_eq_second_order = info_source_data['pi_eq_second_order'];
    const t99_short = info_source_data['t99_short'];
    const t99_long = info_source_data['t99_long'];
    const t99_pl = info_source_data['t99_pl'];
    const t99_pi = info_source_data['t99_pi'];
    const t99_diff = info_source_data['t99_diff'];
    const t99_pseudo = info_source_data['t99_pseudo'];
    const diff_calc_resolution = info_source_data['diff_calc_resolution'];
    const t99Ntimes = info_source_data['t99Ntimes'];

    if (cb_obj.name == 'pro_spinner') {p0[0] = cb_obj.value;}
    if (cb_obj.name == 'lig_spinner') {l0[0] = cb_obj.value;}
    if (cb_obj.name == 'inhib_spinner') {i0[0] = cb_obj.value;}
    if (cb_obj.name == 'lig_koff_spinner') {pl_koff[0] = cb_obj.value;}
    if (cb_obj.name == 'lig_kon_spinner') {pl_kon[0] = cb_obj.value;}
    if (cb_obj.name == 'inhib_koff_spinner') {pi_koff[0] = cb_obj.value;}
    if (cb_obj.name == 'inhib_kon_spinner') {pi_kon[0] = cb_obj.value;}
    const pl_kd_js = pl_koff[0] / pl_kon[0] * 1e9;
    const pi_kd_js = pi_koff[0] / pi_kon[0] * 1e9;
    pl_kd[0] = pl_kd_js;
    pi_kd[0] = pi_kd_js;

    function makeArr(startValue, stopValue, cardinality) {
        var arr = [];
        var step = (stopValue - startValue) / (cardinality - 1);
        for (var i = 0; i < cardinality; i++) {
            arr.push(startValue + (step * i));
            }
        return [arr, step];}

    function second_order_kinetics(p0, l0, koff, kon) {
        const time_list = [];
        const asso_list = [];
        kon = kon * 1e-9;
        const a = 1;
        const b = 0 - (p0 + l0 + koff/kon);
        const c = p0 * l0;
        const pl_1 = (0 - b - Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        const pl_2 = (0 - b + Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        const A = 1 / (pl_1 - pl_2);
        const [num_array, dt] = makeArr(0, 0.99, 100);
        for (var i=0; i< 100; i++){
            const asso_tmp = pl_1 * num_array[i];
            const time_tmp = A / kon * Math.log((pl_2 * (pl_1 - asso_tmp)) / (pl_1 * (pl_2 - asso_tmp)));
            asso_list.push(asso_tmp);
            time_list.push(time_tmp); }
        const t99_second_order = time_list[99];
        return [time_list, asso_list, t99_second_order, pl_1];}

    function second_order_kinetics_rerun(p0, l0, koff, kon, t99) {
        const time_list = [];
        const asso_list = [];
        kon = kon * 1e-9;
        const a = 1;
        const b = 0 - (p0 + l0 + koff/kon);
        const c = p0 * l0;
        const pl_1 = (0 - b - Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        const pl_2 = (0 - b + Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        const [num_loop, dt] = makeArr(0, 1, 100);
        for (var i=0; i< 100; i++){
            const time_tmp = t99 * num_loop[i];
            const C = (pl_1 / pl_2) * Math.exp(kon * (pl_1 - pl_2) * time_tmp);
            const asso_tmp = (C * pl_2 - pl_1) / (C - 1);
            asso_list.push(asso_tmp);
            time_list.push(time_tmp); }
        return [time_list, asso_list];}

    var [time_list_pl_js, asso_list_pl_js, t99_pl_js, pl_eq_second_order_js] = second_order_kinetics(p0[0], l0[0], pl_koff[0], pl_kon[0]);
    var [time_list_pi_js, asso_list_pi_js, t99_pi_js, pi_eq_second_order_js] = second_order_kinetics(p0[0], i0[0], pi_koff[0], pi_kon[0]);
    var t99_short_js;
    var t99_long_js;
    if (t99_pl_js >= t99_pi_js) {
        t99_long_js = t99_pl_js;
        t99_short_js = t99_pi_js;
        [time_list_pi_js, asso_list_pi_js] = second_order_kinetics_rerun(p0[0], i0[0], pi_koff[0], pi_kon[0], t99_pl_js);
    } else {
        t99_long_js = t99_pi_js;
        t99_short_js = t99_pl_js;
        [time_list_pl_js, asso_list_pl_js] = second_order_kinetics_rerun(p0[0], l0[0], pl_koff[0], pl_kon[0], t99_pi_js);
    }
    t99_short[0] = t99_short_js;
    t99_long[0] = t99_long_js;
    t99_pl[0] = t99_pl_js;
    t99_pi[0] = t99_pi_js;
    pl_eq_second_order[0] = pl_eq_second_order_js;
    pi_eq_second_order[0] = pi_eq_second_order_js;
    if (pl_eq_second_order_js > pi_eq_second_order_js) {
        yr.start = -0.05 * pl_eq_second_order_js;
        yr.end = 1.5 * pl_eq_second_order_js;
    } else {
        yr.start = -0.05 * pi_eq_second_order_js;
        yr.end = 1.5 * pi_eq_second_order_js;
        }

    for (var i=0; i< time_list_pl_js.length; i++){
        time_list_pl[i] = time_list_pl_js[i];
        asso_list_pl[i] = asso_list_pl_js[i];
    }
    for (var i=0; i< time_list_pi_js.length; i++){
        time_list_pi[i] = time_list_pi_js[i];
        asso_list_pi[i] = asso_list_pi_js[i];
    }

    function competitive_binding_eq(P0, A0, B0, KA, KB) {
        // KA,  KB units is nM.
        const a = KA + KB + A0 + B0 - P0;
        const b = KB * (A0 - P0) + KA * (B0 - P0) + KA * KB;
        const c = 0 - KA * KB * P0;
        const theta = Math.acos((-2*a**3 + 9*a*b - 27*c) / (2*Math.sqrt((a*a - 3*b)**3)));
        const p_eq = - a / 3 + (2/3) * Math.sqrt(a*a-3*b) * Math.cos(theta/3);
        const pa_eq = p_eq * A0 / (KA + p_eq);
        const pb_eq = p_eq * B0 / (KB + p_eq);
        return [p_eq, pa_eq, pb_eq];}

    function diff_calc_association(p0, l0, i0, pl_koff, pl_kon, pi_koff, pi_kon, pl_comp_eq, t99_short, t99_long, resolution=diff_calc_resolution[0], times=t99Ntimes[0], plot_points=300) {
        pl_kon = pl_kon * 1e-9;
        pi_kon = pi_kon * 1e-9;
        const t99_2 = t99_long * times;
        const [time_loop1, dt1] = makeArr(0, t99_short, resolution);
        const [time_loop2, dt2] = makeArr(t99_short, t99_2, resolution);
        const con_pl_0 = p0 * l0 * pl_kon * dt1;
        const coff_pl_0 = con_pl_0 * pl_koff * dt1;
        const con_pi_0 = p0 * i0 * pi_kon * dt1;
        const coff_pi_0 = con_pi_0 * pi_koff * dt1;
        const con0 = con_pl_0 + con_pi_0;
        const coff0 = coff_pl_0 + coff_pi_0;
        const conlist_p = [con0];
        const cofflist_p = [coff0];
        const consumlist_p = [con0];
        const coffsumlist_p =  [coff0];
        const conlist_l = [con_pl_0];
        const cofflist_l = [coff_pl_0];
        const consumlist_l = [con_pl_0];
        const coffsumlist_l = [coff_pl_0];
        const conlist_i  = [con_pi_0];
        const cofflist_i = [coff_pi_0];
        const consumlist_i  = [con_pi_0];
        const coffsumlist_i = [coff_pi_0];
        for (var i = 2; i < resolution; i++) {
            var con_l = (p0 - consumlist_p[consumlist_p.length-1] + coffsumlist_p[coffsumlist_p.length-1]) * (l0 - consumlist_l[consumlist_l.length-1] + coffsumlist_l[coffsumlist_l.length-1]) * pl_kon * dt1;
            var con_i = (p0 - consumlist_p[consumlist_p.length-1] + coffsumlist_p[coffsumlist_p.length-1]) * (i0 - consumlist_i[consumlist_i.length-1] + coffsumlist_i[coffsumlist_i.length-1]) * pi_kon * dt1;
            var con_p = con_l + con_i;
            conlist_p.push(con_p);
            conlist_l.push(con_l);
            conlist_i.push(con_i);
            consumlist_p.push(consumlist_p[consumlist_p.length-1] + con_p);
            consumlist_l.push(consumlist_l[consumlist_l.length-1] + con_l);
            consumlist_i.push(consumlist_i[consumlist_i.length-1] + con_i);
            var coff_l = (consumlist_l[consumlist_l.length-1] - coffsumlist_l[coffsumlist_l.length-1]) * pl_koff * dt1;
            var coff_i = (consumlist_i[consumlist_i.length-1] - coffsumlist_i[coffsumlist_i.length-1]) * pi_koff * dt1;
            var coff_p = coff_l + coff_i;
            cofflist_p.push(coff_p);
            cofflist_l.push(coff_l);
            cofflist_i.push(coff_i);
            coffsumlist_p.push(coffsumlist_p[coffsumlist_p.length-1] + coff_p);
            coffsumlist_l.push(coffsumlist_l[coffsumlist_l.length-1] + coff_l);
            coffsumlist_i.push(coffsumlist_i[coffsumlist_i.length-1] + coff_i);
        }
        for (var i = 1; i < resolution; i++) {
            var con_l = (p0 - consumlist_p[consumlist_p.length-1] + coffsumlist_p[coffsumlist_p.length-1]) * (l0 - consumlist_l[consumlist_l.length-1] + coffsumlist_l[coffsumlist_l.length-1]) * pl_kon * dt2;
            var con_i = (p0 - consumlist_p[consumlist_p.length-1] + coffsumlist_p[coffsumlist_p.length-1]) * (i0 - consumlist_i[consumlist_i.length-1] + coffsumlist_i[coffsumlist_i.length-1]) * pi_kon * dt2;
            var con_p = con_l + con_i;
            conlist_p.push(con_p);
            conlist_l.push(con_l);
            conlist_i.push(con_i);
            consumlist_p.push(consumlist_p[consumlist_p.length-1] + con_p);
            consumlist_l.push(consumlist_l[consumlist_l.length-1] + con_l);
            consumlist_i.push(consumlist_i[consumlist_i.length-1] + con_i);
            var coff_l = (consumlist_l[consumlist_l.length-1] - coffsumlist_l[coffsumlist_l.length-1]) * pl_koff * dt2;
            var coff_i = (consumlist_i[consumlist_i.length-1] - coffsumlist_i[coffsumlist_i.length-1]) * pi_koff * dt2;
            var coff_p = coff_l + coff_i;
            cofflist_p.push(coff_p);
            cofflist_l.push(coff_l);
            cofflist_i.push(coff_i);
            coffsumlist_p.push(coffsumlist_p[coffsumlist_p.length-1] + coff_p);
            coffsumlist_l.push(coffsumlist_l[coffsumlist_l.length-1] + coff_l);
            coffsumlist_i.push(coffsumlist_i[coffsumlist_i.length-1] + coff_i);
        }
        const pl_list = [];
        pl_list.push(0);
        const time_loop = [];
        for (var i=0; i< resolution; i++){
            time_loop.push(time_loop1[i]);
        }
        for (var i=1; i< resolution; i++){
            time_loop.push(time_loop2[i]);
        }
        const pl_99 = 0.99 * pl_comp_eq;
        const pl_101 = 1.01 * pl_comp_eq;
        var idx_99 = 0;
        for (var i=0; i< 2 * resolution - 2; i++){
            const pl_tmp = consumlist_l[i] - coffsumlist_l[i];
            pl_list.push(pl_tmp);
        }
        for (var i = 2 * resolution - 1; i >= 0; i--) {
            if (pl_list[i] < pl_99 || pl_list[i] > pl_101) {
                idx_99 = i;
                break;
            }
        }

       var t99_approx_idx = 0;
       var diff = Math.abs(t99_long - time_loop[0]);
       for (var i = 0; i < time_loop.length; i++) {
            const newdiff = Math.abs(t99_long - time_loop[i]);
            if (newdiff < diff) {
                diff = newdiff;
                t99_approx_idx = i;
                };
            };

        const [samples_idx_float, dt3] = makeArr(0, t99_approx_idx, plot_points);
        const samples_time = [];
        const samples_pl = [];
        for (var i=0; i< plot_points; i++){
            const samples_idx = Math.floor(samples_idx_float[i]);
            samples_time.push(time_loop[samples_idx]);
            samples_pl.push(pl_list[samples_idx]);
        }
        return [samples_time, samples_pl, time_loop[idx_99], pl_list[pl_list-1]];}

    function kinetics_pseudo_asso(N, L0, I0, k1, k2, k3, k4, t99, plot_points=300) {
        k1 = k1 * 1e-9;
        k3 = k3 * 1e-9;
        const K_A = k1 * L0 + k2;
        const K_B = k3 * I0 + k4;
        const tmp_sqrt = Math.sqrt((K_A - K_B)*(K_A - K_B) + 4 * k1 * k3 * L0 * I0);
        const K_F = 0.5 * (K_A + K_B + tmp_sqrt);
        const K_S = 0.5 * (K_A + K_B - tmp_sqrt);
        const PL_eq = N * k1 * k4 * L0 / (K_F * K_S);
        const PL_99 = 0.99 * PL_eq;
        const PL_101 = 1.01 * PL_eq;
        const t_N = 2;
        const t99_2 = t99 * t_N;
        const [comp_tlist, dt] = makeArr(0, t99_2, plot_points * t_N);
        const PL = [];
        for (var i=0; i< plot_points * t_N; i++){
            const PL_tmp = (N * k1 * L0 / (K_F - K_S)) * (k4 * (K_F - K_S) / (K_F * K_S) + ((k4 - K_F) / K_F) * Math.exp(0 - K_F * comp_tlist[i]) - ((k4 - K_S) / K_S) * Math.exp(0 - K_S * comp_tlist[i]));
            PL.push(PL_tmp);
        }
        var idx_99 = 0;
        for (var i = PL.length - 1; i >= 0; i--) {
            if (PL[i] < PL_99 || PL[i] > PL_101) {
                idx_99 = i;
                break;
            }
        }
        const comp_tlist_out = [];
        const PL_out = [];
        for (var i=0; i< plot_points; i++){
            comp_tlist_out.push(comp_tlist[i]);
            PL_out.push(PL[i]);
        }

        return [comp_tlist_out, PL_out, comp_tlist[idx_99], PL_eq];}

    function ic50_curve(p0, l0, pl_kd, pi_kd) {
        // pl_kd and pi_kd units is nM.
        const a = 1;
        const b = 0 - (p0 + l0 + pl_kd);
        const c = p0 * l0;
        const pl_1 = (0 - b - Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        const [plot_ic50_i_log_list, dt] = makeArr(-4, 6, 100);
        const plot_ic50_i_list = [];
        const plot_ic50_i_free_list = [];
        const plot_ic50_pl_list = [];
        var p_eq, pl_comp_eq, pi_eq;
        for (var i = 0; i < plot_ic50_i_log_list.length; i++) {
            plot_ic50_i_list.push(10**plot_ic50_i_log_list[i]);
            }
        for (var i = 0; i < plot_ic50_i_list.length; i++) {
            [p_eq, pl_comp_eq, pi_eq] = competitive_binding_eq(p0, l0, plot_ic50_i_list[i], pl_kd, pi_kd);
            plot_ic50_pl_list.push(pl_comp_eq);
            plot_ic50_i_free_list.push(plot_ic50_i_list[i] - pi_eq);
            }
        const ic50_approx = pi_kd * (1 + l0 / pl_kd);
        if (ic50_approx > 1e6) {
            return [[], [], [], [], ic50_approx]
            }
        const [inner_ic50_i_list, dt2] = makeArr(0, 100 * ic50_approx, 90000);
        const inner_ic50_pl_list = [];
        const inner_ic50_i_free_list = [];
        for (var i = 0; i < inner_ic50_i_list.length; i++) {
            [p_eq, pl_comp_eq, pi_eq] = competitive_binding_eq(p0, l0, inner_ic50_i_list[i], pl_kd, pi_kd);
            inner_ic50_pl_list.push(pl_comp_eq);
            inner_ic50_i_free_list.push(inner_ic50_i_list[i] - pi_eq);
            }
       var diff = Math.abs(pl_1/2 - inner_ic50_pl_list[0]);
       var ic50_approx_idx = 0;
       for (var i = 0; i < inner_ic50_pl_list.length; i++) {
            const newdiff = Math.abs(pl_1/2 - inner_ic50_pl_list[i]);
            if (newdiff < diff) {
                diff = newdiff;
                ic50_approx_idx = i;
                };
            };
        const apparent_ic50 = inner_ic50_i_list[ic50_approx_idx];
        const ic50 = inner_ic50_i_free_list[ic50_approx_idx];
        return [plot_ic50_i_list, plot_ic50_i_free_list, plot_ic50_pl_list, apparent_ic50, ic50];}

    const [p_eq_js, pl_comp_eq_js, pi_eq_js] = competitive_binding_eq(p0[0], l0[0], i0[0], pl_kd_js, pi_kd_js);
    const [samples_time_js, samples_pl_js, t99_diff_js, pl_end_js] = diff_calc_association(p0[0], l0[0], i0[0], pl_koff[0], pl_kon[0], pi_koff[0], pi_kon[0], pl_comp_eq_js, t99_short_js, t99_long_js);
    const [comp_tlist_js, PL_js, t99_pseudo_js, pl_comp_eq_pseudo_js] = kinetics_pseudo_asso(p0[0], l0[0], i0[0], pl_kon[0], pl_koff[0], pi_kon[0], pi_koff[0], t99_long_js);
    const [plot_ic50_i_list_js, plot_ic50_i_free_list_js, plot_ic50_pl_list_js, apparent_ic50_js, ic50_js] = ic50_curve(p0[0], l0[0], pl_kd_js, pi_kd_js);

    function calc_ki_app4(pt, lt, it, pl_kd) {
        const a = 1;
        const b = 0 - (pt + lt + pl_kd);
        const c = pt * lt;
        const pl_1 = (0 - b - Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        const P0 = pt - pl_1;
        const PL50 = pl_1 / 2;
        const L50 = lt - PL50;
        const I50 = it - pt + pl_kd * PL50 / L50 + PL50;
        const ki = I50 / (L50 / pl_kd + P0 / pl_kd + 1);
        return ki;}
    const calc_ki_js = calc_ki_app4(p0[0], l0[0], apparent_ic50_js, pl_kd_js);

    t99_diff[0] = t99_diff_js;
    t99_pseudo[0] = t99_pseudo_js;
    pl_comp_eq[0] = pl_comp_eq_js;
    pl_comp_eq_pseudo[0] = pl_comp_eq_pseudo_js;
    ic50[0] = ic50_js;
    for (var i=0; i< samples_time_js.length; i++){
        samples_time[i] = samples_time_js[i];
        samples_pl[i] = samples_pl_js[i];
    }
    for (var i=0; i< comp_tlist_js.length; i++){
        comp_tlist[i] = comp_tlist_js[i];
        PL[i] = PL_js[i];
    }
    for (var i=0; i< plot_ic50_i_list_js.length; i++){
        plot_ic50_i_list[i] = plot_ic50_i_list_js[i];
        plot_ic50_i_free_list[i] = plot_ic50_i_free_list_js[i];
        plot_ic50_pl_list[i] = plot_ic50_pl_list_js[i];
    }
info_div.text = '<b><i>K</i><sub>d</sub></b> for Protein and Ligand (<b><i>K</i><sub>d</sub></b> = <b><i>k</i><sub>off-ligand</sub></b>/<b><i>k</i><sub>on-ligand</sub></b>): ' + pl_kd_js.toFixed(5) + ' nM;<br>'+
'<b><i>K</i><sub>d</sub></b> for Protein and Inhibitor (<b><i>K</i><sub>d</sub></b>, or <b><i>K</i><sub>i</sub></b> = <b><i>k</i><sub>off-inhibitor</sub></b>/<b><i>k</i><sub>on-inhibitor</sub></b>): ' + pi_kd_js.toFixed(5) + ' nM;<br>'+
'Exact Time to Reach 99% Equilibrium (<b><i>t</i><sub>0.99</sub></b>), <b>[PL]<sub>eq</sub></b>, and <b>[PI]<sub>eq</sub></b> in the Second-order Association Process of Only Ligand or Inhibitor:<br>' +
'<b><i>t</i><sub>0.99-ligand</sub></b> = ' + t99_pl_js.toExponential(2) +' s; <b><i>t</i><sub>0.99-inhibitor</sub></b> = ' + t99_pi_js.toExponential(2) +' s;<br>' +
'<b>[PL]<sub>eq</sub></b> = ' + pl_eq_second_order_js.toFixed(5) + ' nM; <b>[PI]<sub>eq</sub></b> = ' + pi_eq_second_order_js.toFixed(5) + ' nM;<br>' +
'Approximate <b><i>t</i><sub>0.99</sub></b> in Competitive Second-order Association Process: ' + t99_diff_js.toExponential(2) +' s;<br>' +
'Approximate <b><i>t</i><sub>0.99</sub></b> in Competitive Pseudo-first-order Association Process: ' + t99_pseudo_js.toExponential(2) +' s;<br>' +
'Exact <b>[PL]<sub>eq</sub></b> in Competitive Second-order Binding Process: ' + pl_comp_eq_js.toFixed(5) + ' nM;<br>' +
'Exact <b>[PL]<sub>eq</sub></b> in Competitive Pseudo-first-order Binding Process: ' + pl_comp_eq_pseudo_js.toFixed(5) + ' nM;<br>';

    if (Math.abs(pl_end_js - pl_comp_eq_js) > 0.001 * pl_comp_eq_js) {
        error_div.text = '<p style="background-color: red;">The difference between the simulated <b>[PL]<sub>eq</sub></b> and theoretical <b>[PL]<sub>eq</sub></b> is greater than 0.1% of the theoretical <b>[PL]<sub>eq</sub></b>!</p>';
        }
        
    if (ic50_js < 1e6) {
ic50_div.text = 'Under the Condition of <b>[P]<sub>0</sub></b> = ' + p0[0].toFixed(5) + ' nM, <b>[L]<sub>0</sub></b> = ' + l0[0].toFixed(5) +
' nM, <b><i>K</i><sub>d</sub></b> = ' + pl_kd_js.toFixed(5) +
' nM, and <b><i>K</i><sub>i</sub></b> = ' + pi_kd_js.toFixed(5) + ' nM, <b>IC<sub>50</sub></b> Approximately Equals ' + ic50_js.toFixed(5) + 
' nM, Apparent <b>IC<sub>50</sub></b> Approximately Equals ' + apparent_ic50_js.toFixed(5) + 
' nM by Wang Equation<sup>1</sup>;<br>' +'<b><i>K</i><sub>i</sub></b> Calculated by Equations<sup>2</sup> Equals '+calc_ki_js.toFixed(5) +' nM.';
    } else {
ic50_div.text = 'Under the Condition of <b>[P]<sub>0</sub></b> = ' + p0[0].toFixed(5) + ' nM, <b>[L]<sub>0</sub></b> = ' + l0[0].toFixed(5) +
' nM, <b><i>K</i><sub>d</sub></b> = ' + pl_kd_js.toFixed(5) +
' nM, and <b><i>K</i><sub>i</sub></b> = ' + pi_kd_js.toFixed(5) + ' nM, <b>IC<sub>50</sub></b> is greater than 10 mM by ' +
'Equation <b>IC<sub>50</sub></b>=<b><i>K</i><sub>i</sub></b>(1+<b>[L]<sub>0</sub></b>/<b><i>K</i><sub>d</sub></b>).';
    }

    second_order_pl_source.change.emit();
    second_order_pi_source.change.emit();
    diff_pl_source.change.emit();
    pseudo_pl_source.change.emit();
    ic50_source.change.emit();
    info_source.change.emit();
    """)
    p4_header = Div(text='<h1 style="display:inline; margin-top: 0">Competitive Binding Kinetics - Association</h1>'
                    + '<span>&nbsp;<a href="/binding-curve-viewer/">[Back to Home]</a></span>', height=45)
    pro_div = Div(text='<p>Initial Concentration of Fixed Component (e.g., Protein, <b>[P]<sub>0</sub></b>, nM)'
                       + '<br>0.001 ~ 50000 nM</p>',
                  width=150, height=70, margin=(0, 5, 0, 5))
    lig_div = Div(text='<p>Ligand Conc. (<b>[L]<sub>0</sub></b>, nM)'
                       + '<br>0.001 ~ 50000 nM</p>',
                  width=150, height=35, margin=(0, 5, 0, 5))
    inhib_div = Div(text='<p>Inhibitor Conc. (<b>[I]<sub>0</sub></b>, nM)'
                         + '<br>0.001 ~ 50000 nM</p>',
                    width=150, height=35, margin=(0, 10, 0, 5))
    pro_spinner = Spinner(low=1e-3, high=5e4, value=p0, step=1, title='', name='pro_spinner', width=150)
    pro_spinner.js_on_change('value', comp_callback)
    lig_spinner = Spinner(low=1e-3, high=5e4, value=l0, step=1, title='', name='lig_spinner', width=150)
    lig_spinner.js_on_change('value', comp_callback)
    inhib_spinner = Spinner(low=1e-3, high=5e4, value=i0, step=1, title='', name='inhib_spinner', width=150)
    inhib_spinner.js_on_change('value', comp_callback)
    lig_koff_div = Div(text='<p>Ligand <b><i>k</i><sub>off-ligand</sub></b>'
                            + '<br>1e-6 ~ 1 s<sup>-1</sup></p>',
                       width=150, height=35, margin=(5, 5, 0, 5))
    lig_koff_spinner = Spinner(low=1e-6, high=1, value=pl_koff, step=1e-6, title='', name='lig_koff_spinner', width=150)
    lig_koff_spinner.js_on_change('value', comp_callback)
    lig_kon_div = Div(text='<p>Ligand <b><i>k</i><sub>on-ligand</sub></b>'
                           + '<br>1e3 ~ 1e8 M<sup>-1</sup> s<sup>-1</sup></p>',
                      width=150, height=35, margin=(5, 5, 0, 5))
    lig_kon_spinner = Spinner(low=1e3, high=1e8, value=pl_kon, step=1e3, title='', name='lig_kon_spinner', width=150)
    lig_kon_spinner.js_on_change('value', comp_callback)
    inhib_koff_div = Div(text='<p>Inhibitor <b><i>k</i><sub>off-inhibitor</sub></b>'
                              + '<br>1e-6 ~ 1 s<sup>-1</sup></p>',
                         width=150, height=35, margin=(5, 5, 0, 5))
    inhib_koff_spinner = Spinner(low=1e-6, high=1, value=pi_koff, step=1e-6, title='', name='inhib_koff_spinner',
                                 width=150)
    inhib_koff_spinner.js_on_change('value', comp_callback)
    inhib_kon_div = Div(text='<p>Inhibitor <b><i>k</i><sub>on-inhibitor</sub></b>'
                             + '<br>1e3 ~ 1e8 M<sup>-1</sup> s<sup>-1</sup></p>',
                        width=150, height=35, margin=(5, 5, 0, 5))
    inhib_kon_spinner = Spinner(low=1e3, high=1e8, value=pi_kon, step=1e3, title='', name='inhib_kon_spinner',
                                width=150)
    inhib_kon_spinner.js_on_change('value', comp_callback)

    return p4_header, pro_div, lig_div, inhib_div, \
           pro_spinner, lig_spinner, inhib_spinner, \
           lig_koff_div, lig_koff_spinner, lig_kon_div, \
           lig_kon_spinner, inhib_koff_div, inhib_koff_spinner, \
           inhib_kon_div, inhib_kon_spinner, info_div, error_div, ic50_div, ref_div, \
           plot_kinetics_comp_asso, ic50_fig


def competitive_binding_dissociation(p0=200, l0=150, i0=200, pl_koff=1e-2, pl_kon=1e5, pi_koff=1e-2, pi_kon=1e5,
                                     diff_calc_resolution=400000, t99Ntimes=10, volume_ratio=1.0):
    r"""The simulation of the competitive binding of dissociation is initiated by mixing the equilibrated
    protein and ligand with the inhibitor in a certain volume ratio. p0, l0, and i0 is respectively the total
    concentration of the protein, ligand, and inhibitor.

    .. math:: :label: eq_5.1

        \mathrm{P} + \mathrm{L} \mathrel{\mathop{\rightleftarrows}^{k_{\mathrm{ligand\textnormal{-}off}}}_{k_{\mathrm{ligand\textnormal{-}on}}}} \mathrm{PL}

    .. math:: :label: eq_5.2

        \mathrm{P} + \mathrm{I} \mathrel{\mathop{\rightleftarrows}^{k_{\mathrm{inhibitor\textnormal{-}off}}}_{k_{\mathrm{inhibitor\textnormal{-}on}}}} \mathrm{PI}

    The competitive binding is described by above equations. To describe the dissociation
    process of the competitive binding, we used the the numerical method to simulate the
    second-order kinetic process of the competitive binding of dissociation. See the Methods of :ref:`the
    paper<bcv_paper>` for detailed expression.

    Parameters
    ----------
    p0 : float
        The total concentration of the protein before mixing in the unit of nM (default value is 200).
    l0 : float
        The total concentration of the ligand before mixing in the unit of nM (default value is 150).
    i0 : float
        The total concentration of the inhibitor before mixing in the unit of nM (default value is 150).
    pl_koff : float
        The value of the ligand's dissociation rate constant (:math:`k_{\mathrm{ligand\textnormal{-}off}}`) in the unit of s\ :sup:`-1`
        (default value is 1e-2).
    pl_kon : float
        The value of the ligand's association rate constant (:math:`k_{\mathrm{ligand\textnormal{-}on}}`) in the unit of M\ :sup:`-1` s\ :sup:`-1`
        (default value is 1e5).
    pi_koff : float
        The value of the inhibitor's dissociation rate constant (:math:`k_{\mathrm{inhibitor\textnormal{-}off}}`) in the unit of s\ :sup:`-1`
        (default value is 1e-2).
    pi_kon : float
        The value of the inhibitor's association rate constant (:math:`k_{\mathrm{inhibitor\textnormal{-}on}}`) in the unit of M\ :sup:`-1` s\ :sup:`-1`
        (default value is 1e5).
    diff_calc_resolution : int
        The resolution of the numerical simulation (default value is 400,000).
    t99Ntimes : int
        The value affects the time length of the simulation (default value is 10).
    volume_ratio : float
        The volume ratio of the mixture of protein and ligand to the inhibitor. By default, the ratio is 1.0, i.e.,
        the volume of the mixture of protein and ligand is the same as the inhibitor. Its range is from 1.0 to 1000.0.

    Returns
    -------
    Tuple
        A tuple of widgets and figures, including ``p5_header``, ``pro_div``, ``volume_ratio_div``,
        ``lig_div``, ``inhib_div``, ``pro_spinner``, ``volume_ratio_spinner``, ``lig_spinner``,
        ``inhib_spinner``, `lig_koff_div``, ``lig_koff_spinner``, ``lig_kon_div``,  ``lig_kon_spinner``,
        ``inhib_koff_div``, ``inhib_koff_spinner``, ``inhib_kon_div``, ``inhib_kon_spinner``,
        ``info_div``, ``error_div``, ``ic50_div``, ``ref_div``,  ``plot_kinetics_comp_disso``,
        ``ic50_fig`` can be used to generate the web application.
    """
    def second_order_kinetics(p0, l0, koff, kon):
        kon = kon * 1e-9  # convert the unit from M^-1 s^-1 to nM^-1 s^-1
        a = 1
        b = 0 - (p0 + l0 + koff/kon)
        c = p0 * l0
        pl_1 = (0 - b - np.sqrt(b * b - 4 * a * c)) / (2 * a)  # pl_1 is [PL]eq, pl_1 < pl_2
        pl_2 = (0 - b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
        A = 1 / (pl_1 - pl_2)
        num_array = np.arange(0, 100, 1) / 100
        asso_list = pl_1 * num_array
        time_list = A / kon * np.log((pl_2 * (pl_1 - asso_list)) / (pl_1 * (pl_2 - asso_list)))
        t99 = time_list[99]
        return time_list, asso_list, t99, pl_1

    def calc_disso_t99(koff):
        disso_t99 = np.log(100) / koff
        return disso_t99

    def competitive_binding_eq(P0, A0, B0, KA, KB):
        # The units of P0, A0, and B0 are nM
        # The units of KA and KB are M
        # Equations from FEBS Letters 360 (1995) 111-114
        # KA = koff1 / kon1
        # KB = koff2 / kon2
        KA = KA * 1e9  # Change unit of M to nM
        KB = KB * 1e9
        a = KA + KB + A0 + B0 - P0
        b = KB * (A0 - P0) + KA * (B0 - P0) + KA * KB
        c = -KA * KB * P0
        theta = np.arccos((-2*a**3 + 9*a*b - 27*c) / (2*np.sqrt((a**2 - 3*b)**3)))
        p_eq = -a / 3 + (2/3) * np.sqrt(a**2-3*b) * np.cos(theta/3)
        pa_eq = p_eq * A0 / (KA + p_eq)
        pb_eq = p_eq * B0 / (KB + p_eq)
        return p_eq, pa_eq, pb_eq

    def diff_calc_dissociation(pl0, p0, l0, i0, pl_koff, pl_kon, pi_koff, pi_kon, pl_comp_eq, t99_short, t99_long,
                               resolution=diff_calc_resolution, t99Ntimes=t99Ntimes, plot_points=300):
        pi0 = 0
        # The units of p0, l0, and i0 are nM
        # The units of pl_koff and pi_koff are s-1
        # The units of pl_kon and pi_kon are M-1 s-1
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
        time_loop = np.append(time_loop1, time_loop2[1:])   # length should be 2 * resolution - 1
        pl_99 = 0.99 * pl_comp_eq
        pl_101 = 1.01 * pl_comp_eq
        idx_99 = 0  # the index of approximate pl at 99% equilibrium
        # if idx_99 doesn't change, the starting [PL] is within [pl_99, pl_101]
        # the length of consumlist_l is 2 * resolution - 2
        for i in range(0, 2 * resolution - 2):
            pl_list.append(pl0 + consumlist_l[i] - coffsumlist_l[i])
        for i in reversed(range(2 * resolution - 1)):
            if pl_list[i] < pl_99 or pl_list[i] > pl_101:
                idx_99 = i
                break
        t99_approx_idx = min(range(len(time_loop)), key=lambda i: abs(time_loop[i] - t99_long))
        samples_idx = np.linspace(0, t99_approx_idx, plot_points).astype(int)
        samples_time = [time_loop[idx] for idx in samples_idx]
        samples_pl = [pl_list[idx] for idx in samples_idx]
        return samples_time, samples_pl, time_loop[idx_99]

    def ic50_curve(p0, l0, pl_kd, pi_kd):
        # The units of p0 and l0 are nM
        # The units of pl_kd and pi_kd are nM
        pl_kd_M, pi_kd_M = pl_kd * 1e-9, pi_kd * 1e-9
        # without inhibitor the PL_eq equals pl_1
        a = 1
        b = 0 - (p0 + l0 + pl_kd)
        c = p0 * l0
        pl_1 = (0 - b - np.sqrt(b * b - 4 * a * c)) / (2 * a)
        ic50_approx = pi_kd * (1 + l0 / pl_kd)
        if ic50_approx > 1e6:
            return None, None, None, None, None
        plot_ic50_i_log_list = np.linspace(-4, np.log10(ic50_approx*100), 100)
        plot_ic50_i_list = []
        plot_ic50_i_free_list = []
        plot_ic50_pl_list = []
        for i_log in plot_ic50_i_log_list:
            plot_ic50_i_list.append(10**i_log)
        for i in plot_ic50_i_list:
            p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(p0, l0, i, pl_kd_M, pi_kd_M)
            plot_ic50_pl_list.append(pl_comp_eq)
            plot_ic50_i_free_list.append(i - pi_eq)
        ic50_approx = pi_kd * (1 + l0 / pl_kd)
        inner_ic50_i_list = np.linspace(0, 100 * ic50_approx, 90000)
        inner_ic50_i_free_list = []  #this is I free
        inner_ic50_pl_list = []
        for i in inner_ic50_i_list:
            p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(p0, l0, i, pl_kd_M, pi_kd_M)
            #pl_quan = p0 * l0 / (kd * (1 + i / ki) + l0)
            inner_ic50_pl_list.append(pl_comp_eq)
            inner_ic50_i_free_list.append(i - pi_eq)
        ic50_approx_idx = min(range(len(inner_ic50_pl_list)), key=lambda i: abs(inner_ic50_pl_list[i] - (pl_1 / 2)))
        apparent_ic50 = inner_ic50_i_list[ic50_approx_idx]
        ic50 = inner_ic50_i_free_list[ic50_approx_idx]
        return plot_ic50_i_list, plot_ic50_i_free_list, plot_ic50_pl_list, apparent_ic50, ic50

    time_list_pl, asso_list_pl, t99_pl, pl_eq_second_order = second_order_kinetics(p0, l0, pl_koff, pl_kon)

    p_total_mix = p0 * volume_ratio / (1 + volume_ratio)
    l_total_mix = l0 * volume_ratio / (1 + volume_ratio)
    i_total_mix = i0 / (1 + volume_ratio)
    p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(P0=p_total_mix, A0=l_total_mix, B0=i_total_mix,
                                                     KA=pl_koff/pl_kon, KB=pi_koff/pi_kon)
    pl_mix = pl_eq_second_order * volume_ratio / (1 + volume_ratio)
    p_mix = (p0 - pl_eq_second_order) * volume_ratio / (1 + volume_ratio)
    l_mix = (l0 - pl_eq_second_order) * volume_ratio / (1 + volume_ratio)
    i_mix = i0 / (1 + volume_ratio)

    time_list_pi, asso_list_pi, t99_pi, pi_eq_second_order = second_order_kinetics(p_mix, i_mix, pi_koff, pi_kon)
    pl_disso_t99 = calc_disso_t99(pl_koff)
    t99_short = None
    t99_long = None
    if pl_disso_t99 >= t99_pi:
        t99_long = pl_disso_t99
        t99_short = t99_pi
    else:
        t99_long = t99_pi
        t99_short = pl_disso_t99

    samples_time_disso, samples_pl_disso, t99_diff_disso = diff_calc_dissociation(pl0=pl_mix, p0=p_mix, l0=l_mix, i0=i_mix,
                               pl_koff=pl_koff, pl_kon=pl_kon, pi_koff=pi_koff, pi_kon=pi_kon,
                               pl_comp_eq=pl_comp_eq, t99_short=t99_short, t99_long=t99_long)
    diff_pl_disso = dict(samples_time_disso=samples_time_disso, samples_pl_disso=samples_pl_disso)
    pl_kd = pl_koff / pl_kon * 1e9
    pi_kd = pi_koff / pi_kon * 1e9
    plot_ic50_i_list, plot_ic50_i_free_list, plot_ic50_pl_list, apparent_ic50, ic50 = ic50_curve(p_total_mix, l_total_mix, pl_kd, pi_kd)

    def calc_ki_app4(pt, lt, it, pl_kd):
        a = 1
        b = 0 - (pt + lt + pl_kd)
        c = pt * lt
        pl_1 = (0 - b - np.sqrt(b * b - 4 * a * c)) / (2 * a)
        # use the solutions in Anal. Biochem. 2004, 332 (2), 261-73.
        # pt - pl_1 is the P0 in original paper.
        # pl_1 is the PL0 in original paper.
        P0 = pt - pl_1
        PL50 = pl_1 / 2
        L50 = lt - PL50
        I50 = it - pt + pl_kd * PL50 / L50 + PL50
        ki = I50 / (L50 / pl_kd + P0 / pl_kd + 1)
        return ki
    calc_ki = calc_ki_app4(p_total_mix, l_total_mix, apparent_ic50, pl_kd)

    ic50_data = dict(plot_ic50_i_list=plot_ic50_i_list, plot_ic50_i_free_list=plot_ic50_i_free_list,
                     plot_ic50_pl_list=plot_ic50_pl_list)
    info_data = dict(p0=[p0], l0=[l0], i0=[i0], pl_kd=[pl_kd], pi_kd=[pi_kd], ic50=[ic50], apparent_ic50=[apparent_ic50],
                     pl_koff=[pl_koff], pl_kon=[pl_kon], pi_koff=[pi_koff], pi_kon=[pi_kon],
                     pl_comp_eq=[pl_comp_eq], pl_eq_second_order=[pl_eq_second_order],
                     pi_eq_second_order=[pi_eq_second_order], t99_short=[t99_short], t99_long=[t99_long],
                     t99_pl=[t99_pl], t99_pi=[t99_pi], t99_diff_disso=[t99_diff_disso],
                     volume_ratio=[volume_ratio],p_total_mix=[p_total_mix],
                     l_total_mix=[l_total_mix], i_total_mix=[i_total_mix],
                     pl_mix=[pl_mix], p_mix=[p_mix], l_mix=[l_mix], i_mix=[i_mix],
                     diff_calc_resolution=[diff_calc_resolution], t99Ntimes=[t99Ntimes])
    diff_pl_disso_source = ColumnDataSource(data=diff_pl_disso)
    ic50_source = ColumnDataSource(data=ic50_data)
    info_source = ColumnDataSource(data=info_data)

    plot_kinetics_comp_disso = figure(width=480, height=575, tools=[], output_backend='svg',
                                      y_range=(0.95 * pl_comp_eq, 1.1 * pl_mix), )
    plot_kinetics_comp_disso.toolbar.logo = None
    plot_kinetics_comp_disso.xaxis.axis_label = "Time (s)"
    plot_kinetics_comp_disso.xaxis.axis_label_text_font_style = 'normal'
    plot_kinetics_comp_disso.yaxis.axis_label_text_font_style = 'normal'
    plot_kinetics_comp_disso.yaxis.axis_label = "Concentration of Protein-Ligand Complex ([PL], nM)"
    plot_kinetics_comp_disso.axis.minor_tick_in = -3
    plot_kinetics_comp_disso.axis.minor_tick_out = 6
    comp_line_diff = plot_kinetics_comp_disso.line('samples_time_disso', 'samples_pl_disso',
                                                  source=diff_pl_disso_source, color='#fc8d62',
                                                  line_width=3, line_alpha=0.8, legend_label='Competition-simulation')
    plot_kinetics_comp_disso.legend.location = 'top_left'

    ic50_fig = figure(width=480, height=452, x_axis_type='log', tools=[], output_backend='svg')
    ic50_fig.toolbar.logo = None
    ic50_fig.xaxis.axis_label = "Concentration of Free/Total Inhibitor (nM)"
    ic50_fig.xaxis.axis_label_text_font_style = 'normal'
    ic50_fig.yaxis.axis_label_text_font_style = 'normal'
    ic50_fig.yaxis.axis_label = "Concentration of Protein-Ligand Complex ([PL], nM)"
    ic50_fig.axis.minor_tick_in = -3
    ic50_fig.axis.minor_tick_out = 6
    apparent_ic50_line = ic50_fig.line('plot_ic50_i_list', 'plot_ic50_pl_list', source=ic50_source,
                                       color='#fc8d62', line_width=2, line_alpha=0.8, legend_label='[I] total')
    ic50_line = ic50_fig.line('plot_ic50_i_free_list', 'plot_ic50_pl_list', source=ic50_source,
                              color='#66c2a5', line_width=2, line_alpha=0.8, legend_label='[I] free')
    ic50_fig.legend.location = 'bottom_left'
    info_div = Div(
        text=
        '<b><i>K</i><sub>d</sub></b> for Protein and Ligand (<b><i>K</i><sub>d</sub></b> = <b><i>k</i><sub>off-ligand</sub></b>/<b><i>k</i><sub>on-ligand</sub></b>): ' +'{:.5f}'.format(pl_kd) + ' nM;<br>'+
        '<b><i>K</i><sub>d</sub></b> for Protein and Inhibitor (<b><i>K</i><sub>d</sub></b>, or <b><i>K</i><sub>i</sub></b> = <b><i>k</i><sub>off-inhibitor</sub></b>/<b><i>k</i><sub>on-inhibitor</sub></b>): ' +'{:.5f}'.format(pi_kd) + ' nM;<br>'+
        'Mix in Volume Ratio of '+ '{:.2f}'.format(volume_ratio) + ' (<b>V<sub>protein and ligand</sub></b> / <b>V<sub>inhibitor</sub></b>);<br>' +
        'Before the Mix, <b>[PL]<sub>eq</sub></b> and <b><i>t</i><sub>0.99</sub></b> are ' + '{:.5f}'.format(pl_eq_second_order) + ' nM and ' + '%.2e' % t99_pl + ' s;<br>' +
        'At the Time of Mix, <b>[PL]</b> and <b>[I]</b> are ' + '{:.5f}'.format(pl_mix) + ' nM and ' + '{:.5f}'.format(i_mix) + ' nM;<br>' +
        'After the Mix, <b>[PL]<sub>eq</sub></b> and <b><i>t</i><sub>0.99</sub></b> are ' + '{:.5f}'.format(pl_comp_eq) + ' nM and ' + '%.2e' % t99_diff_disso + ' s.<br>',
        width=500, margin=(15, 25, 5, 5))
    error_div = Div(text='', width=500)
    ic50_div = Div(
        text=
        'Under the Condition of <b>[P]<sub>0</sub></b> = ' + '{:.5f}'.format(p_total_mix) + ' nM, <b>[L]<sub>0</sub></b> = ' + '{:.5f}'.format(l_total_mix) +
        ' nM, <b><i>K</i><sub>d</sub></b> = ' + '{:.5f}'.format(pl_kd) +
        ' nM, and <b><i>K</i><sub>i</sub></b> = ' + '{:.5f}'.format(
            pi_kd) + ' nM, <b>IC<sub>50</sub></b> Approximately Equals ' + '{:.5f}'.format(ic50) +
        ' nM, Apparent <b>IC<sub>50</sub></b> Approximately Equals ' + '{:.5f}'.format(apparent_ic50) +
        ' nM by Wang Equation<sup>1</sup>;<br>' +'<b><i>K</i><sub>i</sub></b> Calculated by Equations<sup>2</sup> Equals '+'{:.5f} nM.'.format(calc_ki),
        width=440, margin=(5, 5, 13, 35))
    ref_div = Div(
        text='<b>References</b><br>' + '[1] FEBS Lett. 2003, 360 (2), 111-114.<br>[2] Anal. Biochem. 2004, 332 (2), 261-273.'
    )
    comp_callback = CustomJS(args=dict(diff_pl_disso_source=diff_pl_disso_source,
                                       ic50_source=ic50_source, info_source=info_source,
                                       info_div=info_div,
                                       error_div=error_div, ic50_div=ic50_div,
                                       yr=plot_kinetics_comp_disso.y_range,), code="""
    const diff_pl_disso_source_data = diff_pl_disso_source.data;
    const samples_time_disso = diff_pl_disso_source_data['samples_time_disso'];
    const samples_pl_disso = diff_pl_disso_source_data['samples_pl_disso'];
    const ic50_source_data = ic50_source.data;
    const plot_ic50_i_list = ic50_source_data['plot_ic50_i_list'];
    const plot_ic50_i_free_list = ic50_source_data['plot_ic50_i_free_list'];
    const plot_ic50_pl_list = ic50_source_data['plot_ic50_pl_list'];
    
    const info_source_data = info_source.data;
    const p0 = info_source_data['p0'];
    const l0 = info_source_data['l0'];
    const i0 = info_source_data['i0'];
    const pl_kd = info_source_data['pl_kd'];
    const pi_kd = info_source_data['pi_kd'];
    const ic50 = info_source_data['ic50'];
    const pl_koff = info_source_data['pl_koff'];
    const pl_kon = info_source_data['pl_kon'];
    const pi_koff = info_source_data['pi_koff'];
    const pi_kon = info_source_data['pi_kon'];
    const pl_comp_eq = info_source_data['pl_comp_eq'];
    const pl_eq_second_order = info_source_data['pl_eq_second_order'];
    const pi_eq_second_order = info_source_data['pi_eq_second_order'];
    const t99_short = info_source_data['t99_short'];
    const t99_long = info_source_data['t99_long'];
    const t99_pl = info_source_data['t99_pl'];
    const t99_pi = info_source_data['t99_pi'];
    const t99_diff = info_source_data['t99_diff_disso'];
    const volume_ratio = info_source_data['volume_ratio'];
    const p_total_mix = info_source_data['p_total_mix'];
    const l_total_mix = info_source_data['l_total_mix'];
    const i_total_mix = info_source_data['i_total_mix'];
    const pl_mix = info_source_data['pl_mix'];
    const p_mix = info_source_data['p_mix'];
    const l_mix = info_source_data['l_mix'];
    const i_mix = info_source_data['i_mix'];
    const diff_calc_resolution = info_source_data['diff_calc_resolution'];
    const t99Ntimes = info_source_data['t99Ntimes'];
    
    if (cb_obj.name == 'pro_spinner') {p0[0] = cb_obj.value;}
    if (cb_obj.name == 'volume_ratio_spinner') {volume_ratio[0] = cb_obj.value;}
    if (cb_obj.name == 'lig_spinner') {l0[0] = cb_obj.value;}
    if (cb_obj.name == 'inhib_spinner') {i0[0] = cb_obj.value;}
    if (cb_obj.name == 'lig_koff_spinner') {pl_koff[0] = cb_obj.value;}
    if (cb_obj.name == 'lig_kon_spinner') {pl_kon[0] = cb_obj.value;}
    if (cb_obj.name == 'inhib_koff_spinner') {pi_koff[0] = cb_obj.value;}
    if (cb_obj.name == 'inhib_kon_spinner') {pi_kon[0] = cb_obj.value;}
    const pl_kd_js = pl_koff[0] / pl_kon[0] * 1e9;
    const pi_kd_js = pi_koff[0] / pi_kon[0] * 1e9;
    pl_kd[0] = pl_kd_js;
    pi_kd[0] = pi_kd_js;

    function makeArr(startValue, stopValue, cardinality) {
        var arr = [];
        var step = (stopValue - startValue) / (cardinality - 1);
        for (var i = 0; i < cardinality; i++) {
            arr.push(startValue + (step * i));
            }
        return [arr, step];}
        
    function calc_disso_t99(koff) {
        const disso_t99 = Math.log(100) / koff;
        return disso_t99;}
        
    function second_order_kinetics(p0, l0, koff, kon) {
        const time_list = [];
        const asso_list = [];
        kon = kon * 1e-9;
        const a = 1;
        const b = 0 - (p0 + l0 + koff/kon);
        const c = p0 * l0;
        const pl_1 = (0 - b - Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        const pl_2 = (0 - b + Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        const A = 1 / (pl_1 - pl_2);
        const [num_array, dt] = makeArr(0, 0.99, 100);
        for (var i=0; i< 100; i++){
            const asso_tmp = pl_1 * num_array[i];
            const time_tmp = A / kon * Math.log((pl_2 * (pl_1 - asso_tmp)) / (pl_1 * (pl_2 - asso_tmp)));
            asso_list.push(asso_tmp);
            time_list.push(time_tmp); }
        const t99_second_order = time_list[99];
        return [time_list, asso_list, t99_second_order, pl_1];}
    
    var [time_list_pl_js, asso_list_pl_js, t99_pl_js, pl_eq_second_order_js] = second_order_kinetics(p0[0], l0[0], pl_koff[0], pl_kon[0]);
    const pl_disso_t99_js = calc_disso_t99(pl_koff[0]);
    
    t99_pl[0] = t99_pl_js;
    t99_pi[0] = t99_pi_js;
    pl_eq_second_order[0] = pl_eq_second_order_js;
    pi_eq_second_order[0] = pi_eq_second_order_js;
    const p_total_mix_js = p0[0] * volume_ratio[0] / (1 + volume_ratio[0]);
    const l_total_mix_js = l0[0] * volume_ratio[0] / (1 + volume_ratio[0]);
    const i_total_mix_js = i0[0] / (1 + volume_ratio[0]);
    const pl_mix_js = pl_eq_second_order_js * volume_ratio[0] / (1 + volume_ratio[0]);
    const p_mix_js = (p0[0] - pl_eq_second_order_js) * volume_ratio[0] / (1 + volume_ratio[0]);
    const l_mix_js = (l0[0] - pl_eq_second_order_js) * volume_ratio[0] / (1 + volume_ratio[0]);
    const i_mix_js = i0[0] / (1 + volume_ratio[0]);
    p_total_mix[0] = p_total_mix_js;
    l_total_mix[0] = l_total_mix_js;
    i_total_mix[0] = i_total_mix_js;
    pl_mix[0] = pl_mix_js;
    p_mix[0] = p_mix_js;
    l_mix[0] = l_mix_js;
    i_mix[0] = i_mix_js;
    
    var [time_list_pi_js, asso_list_pi_js, t99_pi_js, pi_eq_second_order_js] = second_order_kinetics(p_mix[0], i_mix[0], pi_koff[0], pi_kon[0]);
    var t99_short_js;
    var t99_long_js;
    if (pl_disso_t99_js >= t99_pi_js) {
        t99_long_js = pl_disso_t99_js;
        t99_short_js = t99_pi_js;
    } else {
        t99_long_js = t99_pi_js;
        t99_short_js = pl_disso_t99_js;
    }
    t99_short[0] = t99_short_js;
    t99_long[0] = t99_long_js;
    
    function competitive_binding_eq(P0, A0, B0, KA, KB) {
        // KA,  KB units is nM.
        const a = KA + KB + A0 + B0 - P0;
        const b = KB * (A0 - P0) + KA * (B0 - P0) + KA * KB;
        const c = 0 - KA * KB * P0;
        const theta = Math.acos((-2*a**3 + 9*a*b - 27*c) / (2*Math.sqrt((a*a - 3*b)**3)));
        const p_eq = - a / 3 + (2/3) * Math.sqrt(a*a-3*b) * Math.cos(theta/3);
        const pa_eq = p_eq * A0 / (KA + p_eq);
        const pb_eq = p_eq * B0 / (KB + p_eq);
        return [p_eq, pa_eq, pb_eq];}
        
    function diff_calc_dissociation(pl0, p0, l0, i0, pl_koff, pl_kon, pi_koff, pi_kon, pl_comp_eq, t99_short, t99_long, resolution=diff_calc_resolution[0], times=t99Ntimes[0], plot_points=300) {
        const pi0 = 0;
        pl_kon = pl_kon * 1e-9;
        pi_kon = pi_kon * 1e-9;
        const t99_2 = t99_long * times;
        const [time_loop1, dt1] = makeArr(0, t99_short, resolution);
        const [time_loop2, dt2] = makeArr(t99_short, t99_2, resolution);
        const con_pl_0 = p0 * l0 * pl_kon * dt1;
        const coff_pl_0 = (con_pl_0 + pl0) * pl_koff * dt1;
        const con_pi_0 = p0 * i0 * pi_kon * dt1;
        const coff_pi_0 = (con_pi_0 + pi0) * pi_koff * dt1;
        const con0 = con_pl_0 + con_pi_0;
        const coff0 = coff_pl_0 + coff_pi_0;
        const conlist_p = [con0];
        const cofflist_p = [coff0];
        const consumlist_p = [con0];
        const coffsumlist_p =  [coff0];
        const conlist_l = [con_pl_0];
        const cofflist_l = [coff_pl_0];
        const consumlist_l = [con_pl_0];
        const coffsumlist_l = [coff_pl_0];
        const conlist_i  = [con_pi_0];
        const cofflist_i = [coff_pi_0];
        const consumlist_i  = [con_pi_0];
        const coffsumlist_i = [coff_pi_0];
        for (var i = 2; i < resolution; i++) {
            var con_l = (p0 - consumlist_p[consumlist_p.length-1] + coffsumlist_p[coffsumlist_p.length-1]) * (l0 - consumlist_l[consumlist_l.length-1] + coffsumlist_l[coffsumlist_l.length-1]) * pl_kon * dt1;
            var con_i = (p0 - consumlist_p[consumlist_p.length-1] + coffsumlist_p[coffsumlist_p.length-1]) * (i0 - consumlist_i[consumlist_i.length-1] + coffsumlist_i[coffsumlist_i.length-1]) * pi_kon * dt1;
            var con_p = con_l + con_i;
            conlist_p.push(con_p);
            conlist_l.push(con_l);
            conlist_i.push(con_i);
            consumlist_p.push(consumlist_p[consumlist_p.length-1] + con_p);
            consumlist_l.push(consumlist_l[consumlist_l.length-1] + con_l);
            consumlist_i.push(consumlist_i[consumlist_i.length-1] + con_i);
            var coff_l = (pl0 + consumlist_l[consumlist_l.length-1] - coffsumlist_l[coffsumlist_l.length-1]) * pl_koff * dt1;
            var coff_i = (pi0 + consumlist_i[consumlist_i.length-1] - coffsumlist_i[coffsumlist_i.length-1]) * pi_koff * dt1;
            var coff_p = coff_l + coff_i;
            cofflist_p.push(coff_p);
            cofflist_l.push(coff_l);
            cofflist_i.push(coff_i);
            coffsumlist_p.push(coffsumlist_p[coffsumlist_p.length-1] + coff_p);
            coffsumlist_l.push(coffsumlist_l[coffsumlist_l.length-1] + coff_l);
            coffsumlist_i.push(coffsumlist_i[coffsumlist_i.length-1] + coff_i);
        }
        for (var i = 1; i < resolution; i++) {
            var con_l = (p0 - consumlist_p[consumlist_p.length-1] + coffsumlist_p[coffsumlist_p.length-1]) * (l0 - consumlist_l[consumlist_l.length-1] + coffsumlist_l[coffsumlist_l.length-1]) * pl_kon * dt2;
            var con_i = (p0 - consumlist_p[consumlist_p.length-1] + coffsumlist_p[coffsumlist_p.length-1]) * (i0 - consumlist_i[consumlist_i.length-1] + coffsumlist_i[coffsumlist_i.length-1]) * pi_kon * dt2;
            var con_p = con_l + con_i;
            conlist_p.push(con_p);
            conlist_l.push(con_l);
            conlist_i.push(con_i);
            consumlist_p.push(consumlist_p[consumlist_p.length-1] + con_p);
            consumlist_l.push(consumlist_l[consumlist_l.length-1] + con_l);
            consumlist_i.push(consumlist_i[consumlist_i.length-1] + con_i);
            var coff_l = (pl0 + consumlist_l[consumlist_l.length-1] - coffsumlist_l[coffsumlist_l.length-1]) * pl_koff * dt2;
            var coff_i = (pi0 + consumlist_i[consumlist_i.length-1] - coffsumlist_i[coffsumlist_i.length-1]) * pi_koff * dt2;
            var coff_p = coff_l + coff_i;
            cofflist_p.push(coff_p);
            cofflist_l.push(coff_l);
            cofflist_i.push(coff_i);
            coffsumlist_p.push(coffsumlist_p[coffsumlist_p.length-1] + coff_p);
            coffsumlist_l.push(coffsumlist_l[coffsumlist_l.length-1] + coff_l);
            coffsumlist_i.push(coffsumlist_i[coffsumlist_i.length-1] + coff_i);
        }
        const pl_list = [];
        pl_list.push(pl0);
        const time_loop = [];
        for (var i=0; i< resolution; i++){
            time_loop.push(time_loop1[i]);
        }
        for (var i=1; i< resolution; i++){
            time_loop.push(time_loop2[i]);
        }
        const pl_99 = 0.99 * pl_comp_eq;
        const pl_101 = 1.01 * pl_comp_eq;
        var idx_99 = 0;
        for (var i=0; i< 2 * resolution - 2; i++){
            const pl_tmp = pl0 + consumlist_l[i] - coffsumlist_l[i];
            pl_list.push(pl_tmp);
        }
        for (var i = 2 * resolution - 1; i >= 0; i--) {
            if (pl_list[i] < pl_99 || pl_list[i] > pl_101) {
                idx_99 = i;
                break;
            }
        }
        
       var t99_approx_idx = 0;
       var diff = Math.abs(t99_long - time_loop[0]);
       for (var i = 0; i < time_loop.length; i++) {
            const newdiff = Math.abs(t99_long - time_loop[i]);
            if (newdiff < diff) {
                diff = newdiff;
                t99_approx_idx = i;
                };
            };

        const [samples_idx_float, dt3] = makeArr(0, t99_approx_idx, plot_points);
        const samples_time = [];
        const samples_pl = [];
        for (var i=0; i< plot_points; i++){
            const samples_idx = Math.floor(samples_idx_float[i]);
            samples_time.push(time_loop[samples_idx]);
            samples_pl.push(pl_list[samples_idx]);
        }
        return [samples_time, samples_pl, time_loop[idx_99], pl_list[pl_list-1]];}
        
    function ic50_curve(p0, l0, pl_kd, pi_kd) {
        // pl_kd and pi_kd units is nM.
        const a = 1;
        const b = 0 - (p0 + l0 + pl_kd);
        const c = p0 * l0;
        const pl_1 = (0 - b - Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        const [plot_ic50_i_log_list, dt] = makeArr(-4, 6, 100);
        const plot_ic50_i_list = [];
        const plot_ic50_i_free_list = [];
        const plot_ic50_pl_list = [];
        var p_eq, pl_comp_eq, pi_eq;
        for (var i = 0; i < plot_ic50_i_log_list.length; i++) {
            plot_ic50_i_list.push(10**plot_ic50_i_log_list[i]);
            }
        for (var i = 0; i < plot_ic50_i_list.length; i++) {
            [p_eq, pl_comp_eq, pi_eq] = competitive_binding_eq(p0, l0, plot_ic50_i_list[i], pl_kd, pi_kd);
            plot_ic50_pl_list.push(pl_comp_eq);
            plot_ic50_i_free_list.push(plot_ic50_i_list[i] - pi_eq);
            }
        const ic50_approx = pi_kd * (1 + l0 / pl_kd);
        if (ic50_approx > 1e6) {
            return [[], [], [], [], ic50_approx]
            }
        const [inner_ic50_i_list, dt2] = makeArr(0, 100 * ic50_approx, 90000);
        const inner_ic50_pl_list = [];
        const inner_ic50_i_free_list = [];
        for (var i = 0; i < inner_ic50_i_list.length; i++) {
            [p_eq, pl_comp_eq, pi_eq] = competitive_binding_eq(p0, l0, inner_ic50_i_list[i], pl_kd, pi_kd);
            inner_ic50_pl_list.push(pl_comp_eq);
            inner_ic50_i_free_list.push(inner_ic50_i_list[i] - pi_eq);
            }
       var diff = Math.abs(pl_1/2 - inner_ic50_pl_list[0]);
       var ic50_approx_idx = 0;
       for (var i = 0; i < inner_ic50_pl_list.length; i++) {
            const newdiff = Math.abs(pl_1/2 - inner_ic50_pl_list[i]);
            if (newdiff < diff) {
                diff = newdiff;
                ic50_approx_idx = i;
                };
            };
        const apparent_ic50 = inner_ic50_i_list[ic50_approx_idx];
        const ic50 = inner_ic50_i_free_list[ic50_approx_idx];
        return [plot_ic50_i_list, plot_ic50_i_free_list, plot_ic50_pl_list, apparent_ic50, ic50];}

    const [p_eq_js, pl_comp_eq_js, pi_eq_js] = competitive_binding_eq(p_total_mix_js, l_total_mix_js, i_total_mix_js, pl_kd_js, pi_kd_js);
    const [samples_time_js, samples_pl_js, t99_diff_js, pl_end_js] = diff_calc_dissociation(pl_mix_js, p_mix_js, l_mix_js, i_mix_js, pl_koff[0], pl_kon[0], pi_koff[0], pi_kon[0], pl_comp_eq_js, t99_short_js, t99_long_js);
    yr.start = 0.95 * pl_comp_eq_js;
    yr.end = 1.1 * pl_mix_js;
    const [plot_ic50_i_list_js, plot_ic50_i_free_list_js, plot_ic50_pl_list_js, apparent_ic50_js, ic50_js] = ic50_curve(p_total_mix_js, l_total_mix_js, pl_kd_js, pi_kd_js);
    
    function calc_ki_app4(pt, lt, it, pl_kd) {
        const a = 1;
        const b = 0 - (pt + lt + pl_kd);
        const c = pt * lt;
        const pl_1 = (0 - b - Math.sqrt(b * b - 4 * a * c)) / (2 * a);
        const P0 = pt - pl_1;
        const PL50 = pl_1 / 2;
        const L50 = lt - PL50;
        const I50 = it - pt + pl_kd * PL50 / L50 + PL50;
        const ki = I50 / (L50 / pl_kd + P0 / pl_kd + 1);
        return ki;}
    const calc_ki_js = calc_ki_app4(p_total_mix_js, l_total_mix_js, apparent_ic50_js, pl_kd_js);

    
    t99_diff[0] = t99_diff_js;
    pl_comp_eq[0] = pl_comp_eq_js;
    ic50[0] = ic50_js;
    for (var i=0; i< samples_time_js.length; i++){
        samples_time_disso[i] = samples_time_js[i];
        samples_pl_disso[i] = samples_pl_js[i];
    }
    for (var i=0; i< plot_ic50_i_list_js.length; i++){
        plot_ic50_i_list[i] = plot_ic50_i_list_js[i];
        plot_ic50_i_free_list[i] = plot_ic50_i_free_list_js[i];
        plot_ic50_pl_list[i] = plot_ic50_pl_list_js[i];
    }
info_div.text = '<b><i>K</i><sub>d</sub></b> for Protein and Ligand (<b><i>K</i><sub>d</sub></b> = <b><i>k</i><sub>off-ligand</sub></b>/<b><i>k</i><sub>on-ligand</sub></b>): ' + pl_kd_js.toFixed(5) + ' nM;<br>'+
'<b><i>K</i><sub>d</sub></b> for Protein and Inhibitor (<b><i>K</i><sub>d</sub></b>, or <b><i>K</i><sub>i</sub></b> = <b><i>k</i><sub>off-inhibitor</sub></b>/<b><i>k</i><sub>on-inhibitor</sub></b>): ' + pi_kd_js.toFixed(5) + ' nM;<br>'+
'Mix in Volume Ratio of '+ volume_ratio[0].toFixed(2) + ' (<b>V<sub>protein and ligand</sub></b> / <b>V<sub>inhibitor</sub></b>);<br>' +
'Before the Mix, <b>[PL]<sub>eq</sub></b> and <b><i>t</i><sub>0.99</sub></b> are ' + pl_eq_second_order_js.toFixed(5) + ' nM and ' + t99_pl_js.toExponential(2) + ' s;<br>' +
'At the Time of Mix, <b>[PL]</b> and <b>[I]</b> are ' + pl_mix_js.toFixed(5) + ' nM and ' + i_mix_js.toFixed(5) + ' nM;<br>' +
'After the Mix, <b>[PL]<sub>eq</sub></b> and <b><i>t</i><sub>0.99</sub></b> are ' + pl_comp_eq_js.toFixed(5) + ' nM and ' + t99_diff_js.toExponential(2) + ' s.<br>';

    if (Math.abs(pl_end_js - pl_comp_eq_js) > 0.001 * pl_comp_eq_js) {
        error_div.text = '<p style="background-color: red;">The difference between the simulated <b>[PL]<sub>eq</sub></b> and theoretical <b>[PL]<sub>eq</sub></b> is greater than 0.1% of the theoretical <b>[PL]<sub>eq</sub></b>!</p>';
        }
        
    if (ic50_js < 1e6) {
ic50_div.text = 'Under the Condition of <b>[P]<sub>0</sub></b> = ' + p_total_mix_js.toFixed(5) + ' nM, <b>[L]<sub>0</sub></b> = ' + l_total_mix_js.toFixed(5) +
' nM, <b><i>K</i><sub>d</sub></b> = ' + pl_kd_js.toFixed(5) +
' nM, and <b><i>K</i><sub>i</sub></b> = ' + pi_kd_js.toFixed(5) + ' nM, <b>IC<sub>50</sub></b> Approximately Equals ' + ic50_js.toFixed(5) + 
' nM, Apparent <b>IC<sub>50</sub></b> Approximately Equals ' + apparent_ic50_js.toFixed(5) + 
' nM by Wang Equation<sup>1</sup>;<br>' +'<b><i>K</i><sub>i</sub></b> Calculated by Equations<sup>2</sup> Equals '+calc_ki_js.toFixed(5) +' nM.';
    } else {
ic50_div.text = 'Under the Condition of <b>[P]<sub>0</sub></b> = ' + p_total_mix_js.toFixed(5) + ' nM, <b>[L]<sub>0</sub></b> = ' + l_total_mix_js.toFixed(5) +
' nM, <b><i>K</i><sub>d</sub></b> = ' + pl_kd_js.toFixed(5) +
' nM, and <b><i>K</i><sub>i</sub></b> = ' + pi_kd_js.toFixed(5) + ' nM, <b>IC<sub>50</sub></b> is greater than 10 mM by ' +
'Equation <b>IC<sub>50</sub></b>=<b><i>K</i><sub>i</sub></b>(1+<b>[L]<sub>0</sub></b>/<b><i>K</i><sub>d</sub></b>).';
    }

    diff_pl_disso_source.change.emit();
    ic50_source.change.emit();
    info_source.change.emit();
    """)
    p5_header = Div(text='<h1 style="display:inline; margin-top: 0">Competitive Binding Kinetics - Dissociation</h1>'
                    + '<span>&nbsp;<a href="/binding-curve-viewer/">[Back to Home]</a></span>', height=45)
    pro_div = Div(text='<p>Initial Concentration of Fixed Component (e.g., Protein, <b>[P]<sub>0</sub></b>, nM)'
                       + '<br>0.001 ~ 50000 nM</p>',
                  width=150, height=70, margin=(0, 5, 0, 5))
    volume_ratio_div = Div(text='<p>Volume Ratio of the Protein and Ligand to the Inhibitor'
                       + '<br>1.00 ~ 1000.00</p>',
                  width=150, height=45, margin=(40, 5, 0, 5))
    lig_div = Div(text='<p>Ligand Conc. (<b>[L]<sub>0</sub></b>, nM)'
                       + '<br>0.001 ~ 50000 nM</p>',
                  width=150, height=35, margin=(0, 5, 0, 5))
    inhib_div = Div(text='<p>Inhibitor Conc. (<b>[I]<sub>0</sub></b>, nM)'
                       + '<br>0.001 ~ 50000 nM</p>',
                  width=150, height=35, margin=(0, 10, 0, 5))
    pro_spinner = Spinner(low=1e-3, high=5e4, value=p0, step=1, title='', name='pro_spinner', width=150)
    pro_spinner.js_on_change('value', comp_callback)
    volume_ratio_spinner = Spinner(low=1.0, high=1000.0, value=volume_ratio, step=0.1, title='',
                                   name='volume_ratio_spinner', width=150)
    volume_ratio_spinner.js_on_change('value', comp_callback)
    lig_spinner = Spinner(low=1e-3, high=5e4, value=l0, step=1, title='', name='lig_spinner', width=150)
    lig_spinner.js_on_change('value', comp_callback)
    inhib_spinner = Spinner(low=1e-3, high=5e4, value=i0, step=1, title='', name='inhib_spinner', width=150)
    inhib_spinner.js_on_change('value', comp_callback)
    lig_koff_div = Div(text='<p>Ligand <b><i>k</i><sub>off-ligand</sub></b>'
                            + '<br>1e-6 ~ 1 s<sup>-1</sup></p>',
                       width=150, height=35, margin=(5, 5, 0, 5))
    lig_koff_spinner = Spinner(low=1e-6, high=1, value=pl_koff, step=1e-6, title='', name='lig_koff_spinner', width=150)
    lig_koff_spinner.js_on_change('value', comp_callback)
    lig_kon_div = Div(text='<p>Ligand <b><i>k</i><sub>on-ligand</sub></b>'
                           + '<br>1e3 ~ 1e8 M<sup>-1</sup> s<sup>-1</sup></p>',
                      width=150, height=35, margin=(5, 5, 0, 5))
    lig_kon_spinner = Spinner(low=1e3, high=1e8, value=pl_kon, step=1e3, title='', name='lig_kon_spinner', width=150)
    lig_kon_spinner.js_on_change('value', comp_callback)
    inhib_koff_div = Div(text='<p>Inhibitor <b><i>k</i><sub>off-inhibitor</sub></b>'
                       + '<br>1e-6 ~ 1 s<sup>-1</sup></p>',
                  width=150, height=35, margin=(5, 5, 0, 5))
    inhib_koff_spinner = Spinner(low=1e-6, high=1, value=pi_koff, step=1e-6, title='', name='inhib_koff_spinner', width=150)
    inhib_koff_spinner.js_on_change('value', comp_callback)
    inhib_kon_div = Div(text='<p>Inhibitor <b><i>k</i><sub>on-inhibitor</sub></b>'
                       + '<br>1e3 ~ 1e8 M<sup>-1</sup> s<sup>-1</sup></p>',
                  width=150, height=35, margin=(5, 5, 0, 5))
    inhib_kon_spinner = Spinner(low=1e3, high=1e8, value=pi_kon, step=1e3, title='', name='inhib_kon_spinner', width=150)
    inhib_kon_spinner.js_on_change('value', comp_callback)

    return p5_header, pro_div, volume_ratio_div, lig_div, inhib_div,\
           pro_spinner, volume_ratio_spinner, lig_spinner, inhib_spinner,\
           lig_koff_div, lig_koff_spinner, lig_kon_div, \
           lig_kon_spinner, inhib_koff_div, inhib_koff_spinner, \
           inhib_kon_div, inhib_kon_spinner, info_div, error_div, ic50_div, ref_div, \
           plot_kinetics_comp_disso, ic50_fig


def hts_simulation(l0=10, kd=100, i0=10000.0):
    r"""This simulation is used to theoretically access the sensitivity of the competition-based high-throughput
    screening. The default experimental conditions are as follows. The :math:`K_\mathrm{d}` of the protein
    and labeled ligand is 100 nM. In experimental groups, after mixing, the initial working concentration of
    the inhibitor (:math:`[\mathrm{I}]_0`) is 10 μM. After mixing, the concentration of the total ligand
    (:math:`[\mathrm{L}]_\mathrm{total}`) is 10 nM. After mixing, the equilibrium concentration of the protein-ligand
    complex in blank controls (:math:`[\mathrm{PL}]_0`) ranged from 10% to 90% of the
    :math:`[\mathrm{L}]_\mathrm{total}`.

    Parameters
    ----------
    l0 : float
        The total concentration of the labeled ligand after mixing in the unit of nM (default value is 10).
    kd : float
        The binding affinity of the protein and labeled ligand. The value of the :math:`K_\mathrm{d}` in the unit of nM (default value is 100).
    i0 : float
        The working concentration of the inhibitor after mixing in the unit of nM (default value is 10000).

    Returns
    -------
    Tuple
        A tuple of widgets and figures, including ``p6_header``, ``info_table``, ``div_kd``, ``spinner_kd``,
        ``div_inh``, ``spinner_inh``, ``div_lig``, ``spinner_lig``, ``plot_rect`` can be used to generate the web
        application.
    """
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

    ki_Y, pl0_X = np.mgrid[100:10100:100, 0.1 * l0:0.91 * l0:0.01 * l0]
    # starting point (pl) is from 1.0 to 9.0 nM, step size = 0.1
    # In FP exp., pt should be the total concentration of ligand ([L]0) which equals 10 nM.
    # The following equation is from the function `dissociation_constant` in binding_curve_viewer.py
    # and is used to calculate the total concentration of ligand to set the configuration of pl.
    # lt = (pt + kd - pl) * pl / (pt - pl), substitute the pt with 10 (nM), pl with pl0_X
    p0 = (l0 + kd - pl0_X) * pl0_X / (l0 - pl0_X)
    i0_uM = i0/1000
    p_eq, pl_comp_eq, pi_eq = competitive_binding_eq(p0, l0, i0, kd, ki_Y)

    X = pl0_X.ravel()
    Y = ki_Y.ravel()
    Z = pl_comp_eq / pl0_X
    Z = Z.ravel()
    p0 = p0.ravel()
    rect_source = ColumnDataSource(dict(X=X, Y=Y, Z=Z, p0=p0))
    para_source = ColumnDataSource(dict(kd=[kd], i0=[i0], l0=[l0]))

    colors = ["#133a62", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d",
              "#b2182b"]

    plot_rect = figure(width=650, height=550, x_range=(0.1 * l0, 0.9 * l0), y_range=(1, 10000),
                       tools=[HoverTool()],
                       tooltips=[('[PL]0', "@X nM"), ('Ki', "@Y nM"),
                                 ('[PL]eq/[PL]0', "@Z"), ('[P]total', "@p0 nM")],
                       output_backend='svg')

    r = plot_rect.rect(x="X", y="Y", width=0.55 * l0,
                       width_units='screen',
                       height=99, source=rect_source,
                       fill_color=linear_cmap("Z", colors, low=0, high=1),
                       line_color=None, )

    plot_rect.add_layout(r.construct_color_bar(
        major_label_text_font_size="10px",
        ticker=BasicTicker(desired_num_ticks=11),
        formatter=PrintfTickFormatter(format="%.1f"),
        label_standoff=6,
        border_line_color=None,
        padding=5,
        title=r"$$\text{[PL]}_\text{eq}/\text{[PL]}_0$$"
    ), 'right')

    plot_rect.grid.visible = False
    plot_rect.xaxis.axis_label = r"$$\text{[PL]}_0\ (\text{nM})$$"
    plot_rect.yaxis.axis_label = r"$$K_{\text{i}}\ (\text{nM})$$"
    plot_rect.xaxis.axis_label_text_font_style = 'normal'
    plot_rect.yaxis.axis_label_text_font_style = 'normal'
    plot_rect.axis.minor_tick_in = -3
    plot_rect.axis.minor_tick_out = 6
    plot_rect.toolbar.logo = None

    hts_callback = CustomJS(
        args=dict(rect_source=rect_source, para_source=para_source,
                  xr=plot_rect.x_range),
        code="""
        function competitive_binding_eq(P0, A0, B0, KA, KB) {
            // KA,  KB units is nM.
            const a = KA + KB + A0 + B0 - P0;
            const b = KB * (A0 - P0) + KA * (B0 - P0) + KA * KB;
            const c = 0 - KA * KB * P0;
            const theta = Math.acos((-2*a**3 + 9*a*b - 27*c) / (2*Math.sqrt((a*a - 3*b)**3)));
            const p_eq = - a / 3 + (2/3) * Math.sqrt(a*a-3*b) * Math.cos(theta/3);
            const pa_eq = p_eq * A0 / (KA + p_eq);
            const pb_eq = p_eq * B0 / (KB + p_eq);
            return [p_eq, pa_eq, pb_eq];}

        function range(start, stop, step) {
            if (typeof stop == 'undefined') {
                // one param defined
                stop = start;
                start = 0;
            }

            if (typeof step == 'undefined') {
                step = 1;
            }

            if ((step > 0 && start >= stop) || (step < 0 && start <= stop)) {
                return [];
            }

            var result = [];
            for (var i = start; step > 0 ? i < stop : i > stop; i += step) {
                result.push(i);
            }

            return result;
        };

        const rect_source_data = rect_source.data;
        const X = rect_source_data['X'];
        const Y = rect_source_data['Y'];
        const Z = rect_source_data['Z'];
        const p0 = rect_source_data['p0'];

        const para_source_data = para_source.data;
        const kd = para_source_data['kd'];
        const i0 = para_source_data['i0'];
        const l0 = para_source_data['l0'];

        const y_N = 100;
        const x_N = 81;
        if (cb_obj.name == 'spinner_kd') {kd[0] = cb_obj.value;}
        if (cb_obj.name == 'spinner_inh') {i0[0] = cb_obj.value*1000;}
        if (cb_obj.name == 'spinner_lig') {l0[0] = cb_obj.value;}
        xr.start = 0.1 * l0[0];
        xr.end = 0.9 * l0[0];
        const pl0_js = range(0.1*l0[0], 0.91*l0[0], 0.01*l0[0]);
        for (var i=0; i< X.length; i++){
            X[i] = pl0_js[i%x_N];
            }

        var p0_tmp, p_eq_js, pl_comp_eq_js, pi_eq_js; 
        for (var i=0; i< X.length; i++){
            p0_tmp = (l0[0] + kd[0] - X[i]) * X[i] / (l0[0] - X[i]);
            [p_eq_js, pl_comp_eq_js, pi_eq_js] = competitive_binding_eq(p0_tmp, l0[0], i0[0], kd[0], Y[i]);
            Z[i] = pl_comp_eq_js/pl0_js[i%x_N];
            p0[i] = p0_tmp;
        }

        rect_source.change.emit();
        para_source.change.emit();
        """)

    p6_header = Div(
        text='<h1 style="display:inline; margin-top: 0">[PL]<sub>eq</sub>/[PL]<sub>0</sub> in the Primary Screen</h1>'
             + '<span>&nbsp;<a href="/binding-curve-viewer/">[Back to Home]</a></span>', height=40,
        width=800, height_policy='auto', )
    div_kd = Div(text="<p><b><i>K</i><sub>d</sub></b> of Protein and Ligand (nM):" +
                      "<br>0.001 ~ 1000 nM</p>", width=300, height=40, margin=(0, 0, 0, 5))
    spinner_kd = Spinner(low=0.001, high=1000, value=kd, step=1, title='', width=240, name='spinner_kd',
                         margin=(0, 0, 0, 5))
    spinner_kd.js_on_change('value', hts_callback)

    div_lig = Div(text="<p>Concentration of Total Labeled Ligand (<b>[L]<sub>total</sub></b>, nM):" +
                       "<br>1.0 ~ 100.0 nM</p>", width=300, height=40, margin=(10, 0, 0, 5))
    spinner_lig = Spinner(low=1, high=100, value=l0, step=1, title='', width=240,
                          margin=(0, 0, 5, 5), name='spinner_lig')
    spinner_lig.js_on_change('value', hts_callback)

    div_inh = Div(text="<p>Initial Working Concentration of Inhibitor (<b>[I]<sub>0</sub></b>, &mu;M):" +
                       "<br>0.1 ~ 50 &mu;M</p>", width=300, height=40, margin=(10, 0, 0, 5))
    spinner_inh = Spinner(low=0.1, high=50, value=i0_uM, step=0.1, title='', width=240,
                          margin=(0, 0, 5, 5), name='spinner_inh')
    spinner_inh.js_on_change('value', hts_callback)

    info_table = Div(text='<p>By default, in experimental groups, after mixing, the initial working concentration of the inhibitor (<b>[I]<sub>0</sub></b>) was '
                          + '10 μM. The concentration of the total ligand (<b>[L]<sub>total</sub></b>) was 10 nM. The '
                          + 'equilibrium concentration of the protein-ligand complex in blank controls '
                          + '(<b>[PL]<sub>0</sub></b>) ranged from 10% to 90% of the <b>[L]<sub>total</sub></b>.</p>'
                          + '<p>Users can obtain the ratio of the equilibrium concentrations of the protein-ligand complex '
                          + 'in experimental groups (<b>[PL]<sub>eq</sub></b>) to <b>[PL]<sub>0</sub></b> under different '
                          + '<b><i>K</i><sub>i</sub></b> and <b>[PL]<sub>0</sub></b> in real time.</p>', width=280, margin=(10, 5, 5, 5))
    return p6_header, info_table, div_kd, spinner_kd, div_inh, spinner_inh, div_lig, spinner_lig, plot_rect

if __name__ == '__main__':
    #01
    output_file('dissociation-constant.html', title='Binding Curve Viewer')
    p1_header, config_col, tab_conc, tab_frac = dissociation_constant()
    layout1 = layout(column(p1_header, row(config_col, tab_conc, tab_frac)))
    show(layout1)

    #02
    output_file('kinetics-equilibria-dissociation.html', title='Binding Curve Viewer')
    p2_header, p2_config_col, plot_kinetics_on, plot_kinetics_disso = kinetics_association_dissociation()
    layout2 = layout(column(p2_header, row(p2_config_col, plot_kinetics_on, plot_kinetics_disso)))
    show(layout2)

    #03
    output_file('iso-affinity-plot.html', title='Binding Curve Viewer')
    p3_header, fig_points, fig_curves, kin_div, kin_table = iso_affinity_plot()
    layout3 = layout(column(p3_header, row(fig_points, fig_curves, column(kin_div, kin_table))))
    show(layout3)

    #04
    output_file('competitive-binding-association.html', title='Binding Curve Viewer')
    p4_header, pro_div, lig_div, inhib_div, \
    pro_spinner, lig_spinner, inhib_spinner, \
    lig_koff_div, lig_koff_spinner, lig_kon_div, \
    lig_kon_spinner, inhib_koff_div, inhib_koff_spinner, \
    inhib_kon_div, inhib_kon_spinner, info_div, error_div, ic50_div, ref_div, \
    plot_kinetics_comp_asso, ic50_fig = competitive_binding_association()
    layout4 = layout(column(p4_header,
        row(
            column(
                   row(
                       column(pro_div, pro_spinner),
                       column(lig_div, lig_spinner, lig_koff_div, lig_koff_spinner, lig_kon_div, lig_kon_spinner),
                       column(inhib_div, inhib_spinner, inhib_koff_div, inhib_koff_spinner, inhib_kon_div, inhib_kon_spinner)
                   ), info_div, ref_div, error_div),
            plot_kinetics_comp_asso, column(ic50_div, ic50_fig))))
    show(layout4)

    #05
    output_file('competitive-binding-dissociation.html', title='Binding Curve Viewer')
    p5_header, pro_div, volume_ratio_div, lig_div, inhib_div, \
    pro_spinner, volume_ratio_spinner, lig_spinner, inhib_spinner, \
    lig_koff_div, lig_koff_spinner, lig_kon_div, \
    lig_kon_spinner, inhib_koff_div, inhib_koff_spinner, \
    inhib_kon_div, inhib_kon_spinner, info_div, error_div, ic50_div, ref_div, \
    plot_kinetics_comp_disso, ic50_fig = competitive_binding_dissociation()
    layout5 = layout(column(p5_header,
        row(
            column(
                   row(
                       column(pro_div, pro_spinner, volume_ratio_div, volume_ratio_spinner),
                       column(lig_div, lig_spinner, lig_koff_div, lig_koff_spinner, lig_kon_div, lig_kon_spinner),
                       column(inhib_div, inhib_spinner, inhib_koff_div, inhib_koff_spinner, inhib_kon_div, inhib_kon_spinner)
                   ), info_div, ref_div, error_div),
            plot_kinetics_comp_disso, column(ic50_div, ic50_fig))))
    show(layout5)

    #06
    output_file('hts-simulation.html', title='Binding Curve Viewer')
    p6_header, info_table, div_kd, spinner_kd, div_inh, spinner_inh, div_lig, spinner_lig, plot_rect = hts_simulation()
    config_col = column(div_kd, spinner_kd, div_lig, spinner_lig, div_inh, spinner_inh, info_table)
    layout6 = layout(column(p6_header, row(config_col, plot_rect)))
    show(layout6)
