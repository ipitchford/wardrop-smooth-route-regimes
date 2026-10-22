#!/usr/bin/env python3
"""Reproduce the two-panel smooth-transition figure and its plotted data.

Requires only the existing NumPy/SciPy numerical solver and Matplotlib.
Floating-point calculations illustrate the exact theorem; they are not proof.
The cost change is evaluated by independently solving the Wardrop equilibrium
with and without access to the connector, then subtracting the resulting
physical-network social costs. All plotted samples are saved in CSV files.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import tempfile

# A writable cache also makes this script usable in a read-only home directory.
os.environ.setdefault('MPLCONFIGDIR',
                      str(Path(tempfile.gettempdir())/'traffic-matplotlib-cache'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, ScalarFormatter
import numpy as np
import scipy

from numerical_demo import A, B, checked, macro, wardrop


PLUM = '#6C3569'
TEAL = '#367E78'
OCHRE = '#AE8238'
CORAL = '#C7654D'
INK = '#292A2D'
GREY = '#67666B'
LIGHT = '#E5E1E6'
K = 84672/78125


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sample_equilibrium(t):
    loads, f, g = wardrop(t)
    physical_cost = float(np.sum(loads*(A+B*loads**4)))
    route_cost = 4*t*(f+g)
    checked(abs(physical_cost-route_cost)<1e-9,
            'Physical edge and route cost calculations disagree')
    return {
        't':float(t),
        'connector_flow':float(loads[4]),
        'direct_1_load':float(loads[0]),
        'direct_2_load':float(loads[1]),
        'shared_1_load':float(loads[2]),
        'shared_2_load':float(loads[3]),
        'direct_bank_delay':float(f),
        'shared_bank_delay':float(g),
        'equilibrium_cost_from_edges':physical_cost,
        'equilibrium_cost_from_routes':float(route_cost)}


def main():
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=project/'figures')
    parser.add_argument('--evidence-dir', type=Path, default=project/'evidence')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.evidence_dir.mkdir(parents=True, exist_ok=True)

    # Include the exact threshold as a sample, avoiding reliance on its rounding.
    times = np.unique(np.r_[np.linspace(.8, 1.3, 251), 1.])
    flow_rows = [sample_equilibrium(float(t)) for t in times]
    checked(all(row['connector_flow']==0 for row in flow_rows if row['t']<=1),
            'Connector is used before the baseline transition')
    checked(all(row['connector_flow']>0 for row in flow_rows if row['t']>1),
            'Connector activation is missing above the baseline transition')

    cost_rows = []
    for delta in np.geomspace(.02, .15, 65):
        row = sample_equilibrium(float(1+delta))
        t = row['t']
        f0, g0 = macro(4*t, 'F'), macro(4*t, 'G')
        # At a Wardrop equilibrium, total edge cost equals total demand times
        # the minimal path delay for each commodity. This avoids using the
        # asymptotic identity to generate the plotted cost difference.
        baseline = 4*t*(f0+g0)
        difference = baseline-row['equilibrium_cost_from_routes']
        reference = K*float(delta)**5
        checked(difference>0, 'Computed cost improvement is not positive')
        row.update({
            'delta':float(delta),
            'no_connector_direct_delay':float(f0),
            'no_connector_shared_delay':float(g0),
            'no_connector_equilibrium_cost':float(baseline),
            'cost_reduction_E0_minus_E':float(difference),
            'fifth_power_reference':float(reference),
            'reduction_over_fifth_power_reference':float(difference/reference)})
        cost_rows.append(row)

    write_csv(args.evidence_dir/'smooth_transition_flow.csv', flow_rows)
    write_csv(args.evidence_dir/'smooth_transition_cost.csv', cost_rows)

    plt.rcParams.update({
        'font.family':'DejaVu Sans',
        'font.size':11,
        'axes.labelsize':12,
        'axes.titlesize':13,
        'axes.titleweight':'semibold',
        'axes.labelcolor':INK,
        'text.color':INK,
        'xtick.color':GREY,
        'ytick.color':GREY,
        'axes.edgecolor':GREY,
        'axes.linewidth':.8,
        'pdf.fonttype':42,
        'ps.fonttype':42,
        'savefig.facecolor':'white',
        'mathtext.fontset':'dejavusans',
    })
    # Exactly 16:9; explicit margins retain that ratio in the saved images.
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 7.2))
    fig.subplots_adjust(left=.077, right=.975, bottom=.205, top=.765, wspace=.29)
    fig.text(.077, .935, 'Smooth route activation in a quartic traffic network',
             fontsize=20, weight='semibold', va='top', color=PLUM)
    fig.text(.077, .876,
             'Five physical edges  ·  Two origin–destination pairs  ·  Each demand is $4t$',
             fontsize=12, color=GREY, va='top')

    for ax in axes:
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y', which='major', color=LIGHT, linewidth=.8)
        ax.set_axisbelow(True)
        ax.tick_params(axis='both', which='major', length=4, pad=6)

    left, right = axes
    left.set_title('(a) Connector activation', loc='left', pad=17)
    left.axvspan(.8, 1., color=TEAL, alpha=.045, linewidth=0)
    left.axvline(1., color=OCHRE, linewidth=1.35, dashes=(2.5,3))
    left.plot([row['t'] for row in flow_rows],
              [row['connector_flow'] for row in flow_rows],
              color=TEAL, linewidth=2.8, solid_capstyle='round')
    left.set_xlim(.8, 1.3)
    left.set_ylim(-.00013, .00465)
    left.set_xlabel('Demand multiplier, $t$', labelpad=10)
    left.set_ylabel('Connector flow, $x(t)$', labelpad=10)
    left.xaxis.set_major_locator(FixedLocator([.8,.9,1.,1.1,1.2,1.3]))
    left.yaxis.set_major_locator(FixedLocator([0,.001,.002,.003,.004]))
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3,-3))
    left.yaxis.set_major_formatter(formatter)
    left.text(.821, .00065, 'Connector unused\n$x(t)=0$',
              color=TEAL, fontsize=11, linespacing=1.5)
    left.text(1.012, .00424, '$t=1$: route activates',
              color=OCHRE, fontsize=10.5)
    left.text(1.027, .00187,
              '$x(t)\\sim \\dfrac{84}{125}(t-1)^3$\nas $t\\downarrow1$',
              fontsize=11, color=TEAL, linespacing=1.7)

    deltas = np.array([row['delta'] for row in cost_rows])
    reductions = np.array([row['cost_reduction_E0_minus_E'] for row in cost_rows])
    references = np.array([row['fifth_power_reference'] for row in cost_rows])
    right.set_title('(b) Cost change after activation', loc='left', pad=17)
    right.loglog(deltas, reductions, color=PLUM, linewidth=2.6,
                 marker='o', markevery=8, markersize=4.3,
                 markeredgecolor='white', markeredgewidth=.6,
                 label='Computed $E_0(t)-E(t)$')
    right.loglog(deltas, references, color=CORAL, linewidth=2.1,
                 dashes=(5,3), label='Leading term $K\\delta^5$')
    right.set_xlim(.0185, .163)
    right.set_ylim(1.5e-9, 1.6e-4)
    right.set_xlabel('Distance past activation, $\\delta=t-1$', labelpad=10)
    right.set_ylabel('Reduction in equilibrium social cost', labelpad=10)
    right.xaxis.set_major_locator(FixedLocator([.02,.03,.05,.1,.15]))
    right.xaxis.set_major_formatter(FuncFormatter(lambda value, _pos:f'{value:g}'))
    right.xaxis.set_minor_formatter(NullFormatter())
    right.legend(loc='lower right', frameon=False, fontsize=10,
                 handlelength=2.7, labelspacing=.8)
    right.annotate('Fifth-power asymptote\n$K=84672/78125$',
                   xy=(.055, K*.055**5), xytext=(.023, 2.1e-5),
                   fontsize=11, color=CORAL, linespacing=1.55,
                   arrowprops={'arrowstyle':'-', 'color':CORAL,
                               'linewidth':.9, 'connectionstyle':'arc3,rad=-.12'})

    fig.text(.077, .096,
             '$E_0$ is the equilibrium cost with the connector unavailable. '
             'The dashed curve is an asymptotic reference as $\\delta\\downarrow0$.',
             fontsize=10.5, color=GREY, va='top')
    fig.text(.077, .053,
             'Numerical illustration; finite-$\\delta$ corrections remain visible. '
             'The exact proof is independent of these calculations.',
             fontsize=10.5, color=GREY, va='top')

    png = args.output_dir/'smooth_transition.png'
    pdf = args.output_dir/'smooth_transition.pdf'
    fig.savefig(png, dpi=200)
    fig.savefig(pdf, metadata={
        'Title':'Smooth route activation in a quartic traffic network',
        'Subject':'Numerical illustration of cubic connector activation and a fifth-order cost change',
        'Keywords':'Wardrop equilibrium; BPR costs; smooth transition; numerical illustration',
        'CreationDate':None,'ModDate':None})
    plt.close(fig)

    metadata = {
        'evidence_type':'Floating-point numerical illustration, not a proof',
        'cost_evaluation':'Independent macro Wardrop solutions with and without connector; no asymptotic formula used for computed costs',
        'dimensions_inches':[12.8,7.2],
        'png_dpi':200,
        'png_pixels':[2560,1440],
        'flow_samples':len(flow_rows),
        'cost_samples':len(cost_rows),
        'flow_range_t':[.8,1.3],
        'cost_range_delta':[.02,.15],
        'physical_edge_order':['f1','f2','g1','g2','h'],
        'edge_free_flow_constants':A.tolist(),
        'edge_quartic_coefficients':B.tolist(),
        'activation_threshold':1.,
        'connector_cubic_coefficient_exact':'84/125',
        'cost_fifth_power_coefficient_exact':'84672/78125',
        'cost_fifth_power_coefficient_float':K,
        'minimum_computed_cost_reduction':float(reductions.min()),
        'endpoint_reduction_over_reference':[
            cost_rows[0]['reduction_over_fifth_power_reference'],
            cost_rows[-1]['reduction_over_fifth_power_reference']],
        'plotted_data_files':['smooth_transition_flow.csv','smooth_transition_cost.csv'],
        'figure_files':['smooth_transition.png','smooth_transition.pdf'],
        'colours':{'plum':PLUM,'teal':TEAL,'ochre':OCHRE,'coral':CORAL},
        'versions':{'numpy':np.__version__,'scipy':scipy.__version__,
                    'matplotlib':matplotlib.__version__}}
    (args.evidence_dir/'smooth_transition_metadata.json').write_text(
        json.dumps(metadata, indent=2)+'\n')
    print(json.dumps({'status':'passed','figures':[str(png),str(pdf)],
                      'flow_samples':len(flow_rows),'cost_samples':len(cost_rows),
                      'minimum_cost_reduction':float(reductions.min()),
                      'endpoint_reduction_over_reference':
                          metadata['endpoint_reduction_over_reference']},indent=2))


if __name__=='__main__':
    main()
