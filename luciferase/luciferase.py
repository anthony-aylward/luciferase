#!/usr/bin/env python3
#===============================================================================
# luciferase.py
#===============================================================================

"""Helper functions and scripts for luciferase reporter data"""




# Imports ======================================================================

import argparse
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from scipy.stats import ttest_ind




# Constants ====================================================================

JSON_EXAMPLES = '''Examples of luciferase reporter data in JSON format:
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
'''




# Functions ====================================================================

def ttest_indicator(a, b):
    pvalue = ttest_ind(a, b).pvalue
    return (
        '***' if pvalue < 0.001
        else '**' if pvalue < 0.01
        else '*' if pvalue < 0.05
        else 'ns'
    )


def luciferase_barplot(
        luc_data: dict,
        output_file_path: str,
        title=''
):
    """Create a barplot from luciferase reporter data

    The input dict should contain either five items or six items. If it
    contains five items, the bars of the resulting plot will have a 2-2-1
    style. If it contains six items, the bars will have a 2-1-2-1 style.

    Parameters
    ----------
    luc_data : dict
        A dictionary containing the luciferase reporter data points
    output_file_path : str
        Path to the output file
    format : str
        Format of the output file [pdf]
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

    luc_data = pd.DataFrame.from_dict(luc_data).transpose()

    if len(luc_data.index) == 5:
        xrange = [.65, 1.35, 2.65, 3.35, 4.6]
        color = ['royalblue', 'skyblue', 'royalblue', 'skyblue', 'lightgrey']
        sig_line_limits = xrange[:4]
        sig_indicators = tuple(
            ttest_indicator(a, b) for a, b in (
                (luc_data.iloc[0, :], luc_data.iloc[1, :]),
                (luc_data.iloc[2, :], luc_data.iloc[3, :])
            )
        )
    elif len(luc_data.index) == 6:
        xrange = [.65, 1.35, 2.05, 3, 3.7, 4.4]
        color = [
            'royalblue',
            'skyblue',
            'lightgrey',
            'seagreen',
            'lightgreen',
            'lightgrey'
        ]
        sig_line_limits = xrange[:2] + xrange[3:5]
        sig_indicators = tuple(
            ttest_indicator(a, b) for a, b in (
                (luc_data.iloc[0, :], luc_data.iloc[1, :]),
                (luc_data.iloc[3, :], luc_data.iloc[4, :])
            )
        )
    elif len(luc_data.index) == 12:
        xrange = [.65, 1.35, 2.05, 3, 3.7, 4.4, 5.35, 6.05, 6.75, 7.7, 8.4, 9.1]
        color = [
            '#F781BF',
            '#FDDAEC',
            'lightgrey',
            '#984EA3',
            '#DECBE4',
            'lightgrey',
            '#FF7F00',
            '#FED9A6',
            'lightgrey',
            '#E41A1C',
            '#FBB4AE',
            'lightgrey'
        ]
        sig_line_limits = xrange[:2] + xrange[3:5] + xrange[6:8] + xrange[9:11]
        sig_indicators = tuple(
            ttest_indicator(a, b) for a, b in (
                (luc_data.iloc[0, :], luc_data.iloc[1, :]),
                (luc_data.iloc[3, :], luc_data.iloc[4, :]),
                (luc_data.iloc[6, :], luc_data.iloc[7, :]),
                (luc_data.iloc[9, :], luc_data.iloc[10, :]),
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
        sig_indicators[0],
        ha='center',
        va='bottom',
        fontsize=24
    )
    ax1.hlines(
        sig_line_height,
        sig_line_limits[2],
        sig_line_limits[3],
        color='black',
        lw=3
    )
    ax1.text(
        (sig_line_limits[2] + sig_line_limits[3]) / 2,
        sig_ind_height,
        sig_indicators[1],
        ha='center',
        va='bottom',
        fontsize=24
    )
    if len(luc_data.index) == 12:
        ax1.hlines(
            sig_line_height,
            sig_line_limits[4],
            sig_line_limits[5],
            color='black',
            lw=3
        )
        ax1.text(
            (sig_line_limits[4] + sig_line_limits[5]) / 2,
            sig_ind_height,
            sig_indicators[2],
            ha='center',
            va='bottom',
            fontsize=24
        )
        ax1.hlines(
            sig_line_height,
            sig_line_limits[6],
            sig_line_limits[7],
            color='black',
            lw=3
        )
        ax1.text(
            (sig_line_limits[6] + sig_line_limits[7]) / 2,
            sig_ind_height,
            sig_indicators[3],
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
        metavar='<path/to/output.{pdf,png}>',
        help='path to the output file'
    )
    parser.add_argument(
        '--title',
        metavar='<"plot title">',
        default='',
        help='title for the barplot'
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    with open(args.data, 'r') as f:
        luc_data = json.load(f)
    luciferase_barplot(luc_data, args.output, title=args.title)
