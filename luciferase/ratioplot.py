#!/usr/bin/env python3
#===============================================================================
# ratioplot.py
#===============================================================================

"""Plot and compare the ref / alt ratios of two luciferase reporter assays"""




# Imports ======================================================================

import argparse
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from estimateratio import estimate_ratio
from luciferase.luciferase import (
    LIGHT_COLOR_PALETTE, load_data, remove_batch_effect
)
from luciferase.ratiotest import luciferase_ratiotest


# Constants ====================================================================

JSON_EXAMPLE = """Example of luciferase reporter data in JSON format:
{
  "Ref, untreated": [33.2, 30.3, 33.3],
  "Alt, untreated": [19.7, 16.2, 18.3],
  "Empty, untreated": [1.0, 1.0, 1.0],
  "Ref, dex": [149.4, 99.7, 124.5],
  "Alt, dex": [44.6, 37.6, 37.7],
  "Empty, dex": [1.1, 1.0, 0.9]
}

The number of entries in the input JSON should be a multiple of 3
"""




# Functions ====================================================================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Plot and compare the ref / alt ratios of two luciferase reporter '
            'assays'
        ),
        epilog=JSON_EXAMPLE,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'data',
        metavar='<path/to/data.{csv,tsv,json,xls,xlsx}>',
        help='path to a file containing luciferase reporter data'
    )
    parser.add_argument(
        'output',
        metavar='<path/to/output.{pdf,png,svg}>',
        help='path to the output file'
    )
    parser.add_argument(
        '--title',
        metavar='<"plot title">',
        default='',
        help='title for the plot'
    )
    parser.add_argument(
        '--conf',
        metavar='<float>',
        type=float,
        default=0.95,
        help='confidence level for confidence intervals [0.95]'
    )
    parser.add_argument(
        '--xlab',
        metavar='<label>',
        nargs='+',
        help='labels for x axis'
    )
    parser.add_argument(
        '--ylab',
        metavar='<label>',
        default='Ratio',
        help='label for y axis [Alt:Ref ratio]'
    )
    parser.add_argument(
        '--colors',
        metavar='<color>',
        nargs='+',
        default=LIGHT_COLOR_PALETTE,
        help='color palette for plotting'
    )
    parser.add_argument(
        '--invert',
        action='store_true',
        help='invert ratios'
    )
    parser.add_argument(
        '--table',
        metavar='<path/to/file.tsv>',
        help='write a table of the data'
    )
    parser.add_argument(
        '--test',
        action='store_true',
        help='perform a permutation test'
    )
    parser.add_argument(
        '--permutations',
        metavar='<int>',
        type=int,
        default=10_000,
        help='number of permutations for statistical test'
    )
    return parser.parse_args()


def sig_indicator(pvalue):
    """Return a significance indicator string for the result of a t-test.

    Parameters
    ----------
    pvalue: float
        the p-value
    
    Returns
    -------
    str
        `***` if p<0.001, `**` if p<0.01, `*` if p<0.05, `ns` otherwise.
    """

    return (
        '***' if pvalue < 0.001
        else '**' if pvalue < 0.01
        else '*' if pvalue < 0.05
        else 'ns'
    )


def luciferase_ratioplot(
    luc_data: dict,
    output_file_path: str,
    title: str = '',
    conf: float = 0.95,
    xlab=None,
    ylab: str = 'Ratio',
    color_palette=LIGHT_COLOR_PALETTE,
    invert=False,
    table=None,
    test=False,
    permutations=10_000
):
    """Plot and compare allelic ratios from luciferase reporter data

    Parameters
    ----------
    luc_data : dict
        A dictionary or pandas.DataFrame containing the luciferase reporter
        data points
    output_file_path : str
        Path to the output file
    title : str
        Title to add to plot
    conf : float
        Confidence level for confidence intervals
    xlab
        list of labels for the x-axis, or None
    ylab : str
        label for the y-axis
    color_palette
        color pallete to use for bars
    invert : bool
        if True, invert ratios
    table : str
        path to write tabular data
        
    
    Examples
    --------
    import luciferase
    luc_data = {
        'Ref, untreated': [33.2, 30.3, 33.3],
        'Alt, untreated': [19.7, 16.2, 18.3],
        'Empty, untreated': [1.0, 1.0, 1.0],
        'Ref, dex': [149.4, 99.7, 124.5],
        'Alt, dex': [44.6, 37.6, 37.7],
        'Empty, dex': [1.1, 1.0, 0.9]
    }
    luciferase.luciferase_ratioplot(
        luc_data,
        'untreated-v-dex.pdf',
        title='DEX v untreated'
    )
    """

    if isinstance(luc_data, dict):
        luc_data = pd.DataFrame.from_dict(luc_data).transpose()
    if 'Batch' in luc_data.index:
        luc_data = remove_batch_effect(luc_data)

    n_groups = int(len(luc_data.index) / 3)
    ratio_data = pd.DataFrame(
        estimate_ratio(
            luc_data.iloc[i + invert],
            luc_data.iloc[i + 1 - invert],
            conf=conf
        )
        for i in range(0, int(len(luc_data.index)), 3)
    )
    if table:
        ratio_data.transpose().to_csv(table, sep='\t')
    ratio_data['xrange'] = [.65 + .7 * x for x in range(n_groups)]
    if not xlab:
        xlab = ['' for _ in range(n_groups)]

    ratio_data['ci_lo'] = [ci[0] for ci in ratio_data['ci']]
    ratio_data['ci_hi'] = [ci[1] for ci in ratio_data['ci']]

    sns.set(font_scale=1.5)
    plt.style.use('seaborn-white')
    fig, ax1 = plt.subplots(1, 1, figsize=(3, 5), dpi=100)
    bars = ax1.bar(
        ratio_data['xrange'], ratio_data['r'], edgecolor='black', lw=2,
        color=color_palette[:n_groups], width=.6
    )
    ax1.vlines(
        ratio_data['xrange'], ratio_data['ci_lo'], ratio_data['ci_hi'],
        color='black', lw=2
    )
    ax1.hlines(
        ratio_data['ci_lo'], ratio_data['xrange'] - 0.1,
        ratio_data['xrange'] + 0.1, color='black', lw=2
    )
    ax1.hlines(
        ratio_data['ci_hi'], ratio_data['xrange'] - 0.1,
        ratio_data['xrange'] + 0.1, color='black', lw=2
    )
    if test:
        pvalue = luciferase_ratiotest(luc_data, permutations)
        sig_line_limits = [x for x in ratio_data['xrange'][:2]]
        max_bar_height = max(ratio_data['ci_hi'])
        sig_line_height = max_bar_height * 1.1
        sig_ind_height = max_bar_height * 1.15
        ax1.hlines(
            sig_line_height,
            sig_line_limits[0],
            sig_line_limits[1],
            color='black',
            lw=3
        )
        ax1.text(
            (sig_line_limits[0] + sig_line_limits[1]) / 2,
            sig_ind_height,
            sig_indicator(pvalue),
            ha='center',
            va='bottom',
            fontsize=24
        )
    ax1.set_xticks(ratio_data['xrange'])
    sns.despine(trim=True, offset=10)
    ax1.tick_params(axis='both', length=6, width=1.25, bottom=True, left=True)
    ax1.set_xticklabels(xlab, rotation=45, ha='right')
    ax1.set_ylabel(ylab, fontsize=20)
    ax1.set_title(title, fontsize=24, y=1.1)
    plt.savefig(output_file_path, bbox_inches='tight')


def main():
    args = parse_arguments()
    luc_data = load_data(args.data)
    luciferase_ratioplot(
        luc_data, args.output, title=args.title, conf=args.conf,
        xlab=args.xlab, ylab=args.ylab, color_palette=args.colors,
        invert=args.invert, table=args.table, test=args.test,
        permutations=args.permutations
    )
