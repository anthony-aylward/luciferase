#!/usr/bin/env python3
#===============================================================================
# luciferase.py
#===============================================================================

"""Helper functions and scripts for luciferase reporter data"""




# Imports ======================================================================

import argparse
import json
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import numpy as np

from itertools import chain
from random import sample
from scipy.stats import ttest_ind
from statistics import mean




# Constants ====================================================================

JSON_EXAMPLES = """Examples of luciferase reporter data in JSON format:
{
  "Non-risk, Fwd": [8.354, 12.725, 8.506],
  "Risk, Fwd": [5.078, 5.038, 5.661],
  "Non-risk, Rev": [9.564, 9.692, 12.622], 
  "Risk, Rev": [10.777, 11.389, 10.598],
  "Empty": [1.042, 0.92, 1.042]
}
{
  "Alt, MIN6": [5.47, 7.17, 6.15],
  "Ref, MIN6": [3.16, 3.04, 4.34],
  "Empty, MIN6": [1.07, 0.83, 0.76],
  "Alt, ALPHA-TC6": [2.50, 3.47, 3.33],
  "Ref, ALPHA-TC6": [2.01, 1.96, 2.31],
  "Empty, ALPHA-TC6": [1.042, 0.92, 1.042]
}

The input JSON should contain either five entries or six entries. If it contains
five entries, the bars of the resulting plot will have a 2-2-1 style. If it
contains six entries, the bars will have a 2-1-2-1 style.

Significance indicators will be written above the bars: `***` if p<0.001,
`**` if p<0.01, `*` if p<0.05, `ns` otherwise.
"""

if os.environ.get('LUCIFERASE_DARK_COLORS'):
    DARK_COLOR_PALETTE = os.environ.get('LUCIFERASE_DARK_COLORS')
else:
    DARK_COLOR_PALETTE = [
        '#F781BF', '#984EA3', '#FF7F00', '#E41A1C', '#377EB8', '#4DAF4A',
        'royalblue', 'seagreen'
    ]

if os.environ.get('LUCIFERASE_LIGHT_COLORS'):
    LIGHT_COLOR_PALETTE = os.environ.get('LUCIFERASE_LIGHT_COLORS')
else:
    LIGHT_COLOR_PALETTE = [
        '#FDDAEC', '#DECBE4', '#FED9A6', '#FBB4AE', '#B3CDE3', '#CCEBC5', 'skyblue', 'lightgreen'
    ]

EMPTY_COLOR = 'lightgrey'




# Functions ====================================================================

def scale_factor_lstsq(x, y):
    """Compute a least-squares best-fit scale factor to transform x to the
    scale of y.

    Parameters
    ----------
    x
        data frame containing independent variable
    y
        data frame containing dependent variable

    Returns
    -------
    float
        The scale factor that best fits x to y
    """

    y_flat = pd.concat([y] * len(x.columns), axis=1).values.flatten()
    x_flat = x.values.flatten()
    return mean(y_i / x_i for x_i, y_i in zip(x_flat, y_flat))


def remove_batch_effect(luc_data):
    """Remove batch effects from luciferase data

    Parameters
    ----------
    luc_data
        data frame containing luciferase data points and batch annotations

    Returns
    -------
    DataFrame
        luciferase data re-normalized to remove batch effects
    """

    if isinstance(luc_data, dict):
        luc_data = pd.DataFrame.from_dict(luc_data).transpose()
    
    construct_indices = [
        i
        for pair in (
            (index, index + 1) for index in range(
                0, int(len(luc_data.index) - 1), 3
            )
        )
        for i in pair
    ]
    construct_data = luc_data.drop('Batch').iloc[construct_indices]
    batch = tuple(int(x) for x in luc_data.loc['Batch',:])
    construct_mean = construct_data.mean(axis=1)
    construct_by_batch = {
        b: construct_data.iloc[
            :, [
                col for col in range(len(construct_data.columns))
                if batch[col] == b
            ]
        ]
        for b in set(batch)
    }
    scale_factors = {
        b: scale_factor_lstsq(x, construct_mean)
        for b, x in construct_by_batch.items()
    }
    scale_factor_mean = mean(scale_factors.values())
    normalized_scale_factor_row = tuple(
        scale_factors[b] / scale_factor_mean for b in batch
    )
    return luc_data.drop('Batch').multiply(normalized_scale_factor_row, axis=1)


def ttest_indicator(a, b, batch=None):
    """Return a significance indicator string for the result of a t-test.

    Parameters
    ----------
    a
        iterable of measurements from population A
    b
        iterable of measurements from population B
    
    Returns
    -------
    str
        `***` if p<0.001, `**` if p<0.01, `*` if p<0.05, `ns` otherwise.
    """

    pvalue = ttest_ind(a, b).pvalue
    print(pvalue)
    return (
        '***' if pvalue < 0.001
        else '**' if pvalue < 0.01
        else '*' if pvalue < 0.05
        else 'ns'
    )


def luciferase_barplot(
    luc_data,
    output_file_path: str,
    title: str = '',
    dark_color_palette=DARK_COLOR_PALETTE,
    light_color_palette=LIGHT_COLOR_PALETTE
):
    """Create a barplot from luciferase reporter data

    The input dict should contain either five items or six items. If it
    contains five items, the bars of the resulting plot will have a 2-2-1
    style. If it contains six items, the bars will have a 2-1-2-1 style.

    Parameters
    ----------
    luc_data
        A dict or pandas.DataFrame containing the luciferase reporter
        data points
    output_file_path : str
        Path to the output file
    title : str
        Title to add to plot
    
    Examples
    --------
    import luciferase
    luc_data = {
        'Non-risk, Fwd': [8.354, 12.725, 8.506],
        'Risk, Fwd': [5.078, 5.038, 5.661],
        'Non-risk, Rev': [9.564, 9.692, 12.622], 
        'Risk, Rev': [10.777, 11.389, 10.598],
        'Empty': [1.042, 0.92, 1.042]
    }
    luciferase.luciferase_barplot(luc_data, 'rs7795896.pdf', title='rs7795896')
    luc_data = {
        'Alt, MIN6': [5.47, 7.17, 6.15],
        'Ref, MIN6': [3.16, 3.04, 4.34],
        'Empty, MIN6': [1.07, 0.83, 0.76],
        'Alt, ALPHA-TC6': [2.50, 3.47, 3.33],
        'Ref, ALPHA-TC6': [2.01, 1.96, 2.31],
        'Empty, ALPHA-TC6': [1.042, 0.92, 1.042]
    }
    luciferase.luciferase_barplot(
        luc_data,
        'min6-v-alpha.pdf',
        title='MIN6 v.Alpha'
    )
    """
    
    if isinstance(luc_data, dict):
        luc_data = pd.DataFrame.from_dict(luc_data).transpose()
    if 'Batch' in luc_data.index:
        luc_data = remove_batch_effect(luc_data)
    
    if 'empty' in luc_data.index[2].casefold():
        n_groups = int(len(luc_data.index) / 3)
        xrange = [
            i * 2.35 + x for i in range(n_groups) for x in (.65, 1.35, 2.05)
        ]
        color = [
            c for i in range(n_groups) for c in (
                dark_color_palette[i], light_color_palette[i], EMPTY_COLOR
            )
        ]
        sig_line_limits = [
            x
            for i in range(0, int(len(luc_data.index)), 3)
            for x in xrange[i:i + 2]
        ]
        sig_indicators = tuple(
            ttest_indicator(a, b) for a, b in (
                (luc_data.iloc[i], luc_data.iloc[i + 1])
                for i in range(0, int(len(luc_data.index)), 3)
            )
        )
    elif 'empty' in luc_data.index[4].casefold():
        n_groups = int(len(luc_data.index) / 5)
        xrange = [
            i * 5.9 + x for i in range(n_groups) for x in (
                .65, 1.35, 2.65, 3.35, 4.6
            )
        ]
        sig_line_limits = [
            x
            for i in range(0, int(len(luc_data.index)), 5)
            for x in xrange[i:i + 4]
        ]
        sig_indicators = tuple(
            ttest_indicator(a, b) for a, b in (
                pair
                for i in range(0, int(len(luc_data.index)), 3)
                for pair in (
                    (luc_data.iloc[i], luc_data.iloc[i + 1]),
                    (luc_data.iloc[i + 2], luc_data.iloc[i + 3]),
                )
            )
        )
    luc_data['mean'] = luc_data.mean(axis=1)
    luc_data['std'] = luc_data.iloc[:,:3].std(axis=1)
    luc_data['xrange'] = xrange

    sns.set(font_scale=1.5)
    plt.style.use('seaborn-white')
    fig, ax1 = plt.subplots(1, 1, figsize=(7, 5), dpi=100)
    bars = ax1.bar(
        luc_data['xrange'],
        luc_data['mean'],
        edgecolor='black',
        lw=2,
        color=color,
        width=.6
    )
    ax1.vlines(
        xrange,
        luc_data['mean'],
        luc_data['mean'] + luc_data['std'],
        color='black',
        lw=2
    )
    ax1.hlines(
        luc_data['mean'] + luc_data['std'],
        luc_data['xrange'] - 0.1,
        luc_data['xrange'] + 0.1,
        color='black',
        lw=2
    )
    
    max_bar_height = max(luc_data['mean'] + luc_data['std'])
    sig_line_height = max_bar_height * 1.1
    sig_ind_height = max_bar_height * 1.15
    for i in range(0, len(sig_line_limits), 2):
        ax1.hlines(
            sig_line_height,
            sig_line_limits[i],
            sig_line_limits[i + 1],
            color='black',
            lw=3
        )
        ax1.text(
            (sig_line_limits[i] + sig_line_limits[i + 1]) / 2,
            sig_ind_height,
            sig_indicators[int(i / 2)],
            ha='center',
            va='bottom',
            fontsize=24
        )

    ax1.set_xticks(xrange)
    sns.despine(trim=True, offset=10)
    ax1.tick_params(axis='both', length=6, width=1.25, bottom=True, left=True)
    ax1.set_xticklabels(list(luc_data.index), rotation=45, ha='right')
    ax1.set_ylabel('F$_{luc}$:R$_{luc}$ ratio', fontsize=20)
    ax1.set_title(title, fontsize=24, y=1.1)

    plt.savefig(output_file_path, bbox_inches='tight')


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Create a barplot from a JSON file containing luciferase reporter'
            ' data'
        ),
        epilog=JSON_EXAMPLES,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'data',
        metavar='<path/to/data.json>',
        help='path to a JSON file containing luciferase reporter data'
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
        help='title for the barplot'
    )
    parser.add_argument(
        '--dark-colors',
        metavar='<color>',
        nargs='+',
        default=DARK_COLOR_PALETTE,
        help='palette of dark colors'
    )
    parser.add_argument(
        '--light-colors',
        metavar='<color>',
        nargs='+',
        default=LIGHT_COLOR_PALETTE,
        help='palette of dark colors'
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    with open(args.data, 'r') as f:
        luc_data = json.load(f)
    luciferase_barplot(
        luc_data,
        args.output,
        title=args.title,
        dark_color_palette=args.dark_colors,
        light_color_palette=args.light_colors
    )
