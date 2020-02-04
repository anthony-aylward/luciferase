#!/usr/bin/env python3
#===============================================================================
# swarmplot.py
#===============================================================================

"""Generate a swarmplot from luciferase reporter data"""

# Imports ======================================================================

import argparse
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from luciferase.luciferase import (
    LIGHT_COLOR_PALETTE, DARK_COLOR_PALETTE, EMPTY_COLOR, JSON_EXAMPLES,
    load_data, remove_batch_effect, ttest_indicator
)



# Functions ====================================================================

def luciferase_swarmplot(
    luc_data,
    output_file_path: str,
    title: str = '',
    dark_color_palette=DARK_COLOR_PALETTE,
    light_color_palette=LIGHT_COLOR_PALETTE,
    table=None,
    transpose=False
):
    """Create a swarmplot from luciferase reporter data

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
    dark_color_palette
        Iterable of dark colors
    light_color_palette
        Iterable of light colors
    table : str
        Filename to write output table
    transpose : bool
        If True, transpose alleles with treatments for significance testing
    
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
    luciferase.luciferase_swarmplot(
        luc_data,
        'rs7795896.pdf',
        title='rs7795896'
    )
    luc_data = {
        'Ref, MIN6': [3.16, 3.04, 4.34],
        'Alt, MIN6': [5.47, 7.17, 6.15],
        'Empty, MIN6': [1.07, 0.83, 0.76],
        'Ref, ALPHA-TC6': [2.01, 1.96, 2.31],
        'Alt, ALPHA-TC6': [2.50, 3.47, 3.33],
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
        xrange = list(range(n_groups * 3))
        if transpose:
            luc_data = luc_data.reindex(
                [
                    luc_data.index[3*j + i]
                    for i in range(3)
                    for j in range(n_groups)
                ]
            )
            color = (
                dark_color_palette[:n_groups]
                + light_color_palette[:n_groups]
                + n_groups * [EMPTY_COLOR]
            )
            sig_line_limits = [
                x
                for i in tuple(range(0, int(len(luc_data.index)), 2))[:-1]
                for x in xrange[i:i + 2]
            ]
            sig_indicators = tuple(
                ttest_indicator(a, b) for a, b in (
                    (luc_data.iloc[i], luc_data.iloc[i + 1])
                    for i in tuple(range(0, int(len(luc_data.index)), 2))[:-1]
                )
            )
        else:
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
        xrange = list(range(n_groups * 5))
        color = [
            c for i in range(n_groups) for c in (
                dark_color_palette[i], light_color_palette[i],
                dark_color_palette[i], light_color_palette[i], EMPTY_COLOR
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
    melted_data = luc_data.transpose().melt()
    if table:
        luc_data.transpose().to_csv(table, sep='\t', index=False)
    m, sd = luc_data.mean(axis=1), luc_data.std(axis=1)
    luc_data['mean'], luc_data['std'] = m, sd
    luc_data['xrange'] = xrange

    sns.set(font_scale=1.5)
    plt.style.use('seaborn-white')
    fig, ax1 = plt.subplots(1, 1, figsize=(7, 5), dpi=100)
    bars = ax1.bar(
        luc_data['xrange'],
        luc_data['mean'],
        edgecolor='black',
        lw=2,
        color='white',
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

    # ax1.set_xticks(xrange)
    sns.swarmplot(
        x='variable',
        y='value',
        data=melted_data,
        palette=color,
        linewidth=1,
        size=8
    )
    sns.despine(trim=True, offset=10)
    ax1.tick_params(axis='both', length=6, width=1.25, bottom=True, left=True)
    ax1.set_xticklabels(list(luc_data.index), rotation=45, ha='right')
    ax1.set_xlabel('')
    ax1.set_ylabel('F$_{luc}$:R$_{luc}$ ratio', fontsize=20)
    ax1.set_title(title, fontsize=24, y=1.1)

    plt.savefig(output_file_path, bbox_inches='tight')


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Create a swarmplot from a JSON file containing luciferase reporter'
            ' data'
        ),
        epilog=JSON_EXAMPLES,
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
    parser.add_argument(
        '--table',
        metavar='<path/to/file.tsv>',
        help='write a table of the data'
    )
    parser.add_argument(
        '--transpose',
        action='store_true',
        help='Transpose alleles with treatments'
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    luc_data = load_data(args.data)
    luciferase_swarmplot(
        luc_data,
        args.output,
        title=args.title,
        dark_color_palette=args.dark_colors,
        light_color_palette=args.light_colors,
        table=args.table,
        transpose=args.transpose
    )
