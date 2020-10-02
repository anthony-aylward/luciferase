#!/usr/bin/env python3
#===============================================================================
# ratiotest.py
#===============================================================================

"""Compare the ref / alt ratios of two luciferase reporter assays with a
permutation test"""




# Imports ======================================================================

import argparse
import itertools
import json
import pandas as pd

from estimateratio import estimate_ratio
from math import log2
from random import sample

from luciferase.luciferase import load_data, remove_batch_effect




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
        '--permutations',
        metavar='<int>',
        type=int,
        default=10_000,
        help='number of permutations'
    )
    return parser.parse_args()


def shuffled_rows(luc_data_segment):
    ncol = len(luc_data_segment.columns)
    pop = sample(
        tuple(
            itertools.chain(
                luc_data_segment.iloc[0, :], luc_data_segment.iloc[1, :]
            )
        ),
        2 * ncol
    )
    empty = tuple(luc_data_segment.iloc[2, :])
    return pd.DataFrame(
        [pop[:ncol], pop[ncol:], empty], index=luc_data_segment.index
    )


def shuffle_luc_data(luc_data):
    return pd.concat(
        shuffled_rows(luc_data.iloc[i:i+3,:]) for i in range(0, int(len(luc_data.index)), 3)
    )


def calculate_ratios(luc_data):
    return pd.DataFrame(
        estimate_ratio(luc_data.iloc[i], luc_data.iloc[i + 1])
        for i in range(0, int(len(luc_data.index)), 3)
    ) 


def test_statistic(ratio_data):
    return abs(log2(ratio_data.loc[0, 'r']) - log2(ratio_data.loc[1, 'r']))


def luciferase_ratiotest(luc_data: dict, permutations: int):
    """Plot and compare allelic ratios from luciferase reporter data

    Parameters
    ----------
    luc_data : dict
        A dictionary or pandas.DataFrame containing the luciferase reporter
        data points
    
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
    luciferase.luciferase_ratioplot(luc_data)
    """

    if isinstance(luc_data, dict):
        luc_data = pd.DataFrame.from_dict(luc_data).transpose()
    if 'Batch' in luc_data.index:
        luc_data = remove_batch_effect(luc_data)
    empirical_dist = tuple(
        test_statistic(calculate_ratios(shuffle_luc_data(luc_data)))
        for _ in range(permutations)
    )
    observed = test_statistic(calculate_ratios(luc_data))
    return sum(x >= observed for x in empirical_dist) / len(empirical_dist)


def main():
    args = parse_arguments()
    luc_data = load_data(args.data)
    pval = luciferase_ratiotest(luc_data, args.permutations)
    print(pval)
