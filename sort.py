#!/usr/bin/env python3
"""Sorts CSVs produced by the experiment."""

import pandas as pd


def sort(fname, outcsv, outdf):

    print('Reading CSV')
    series = pd.read_csv(fname, index_col=[0, 1], squeeze=True,
                         skipinitialspace=True, header=None)

    print('Sorting Series')
    series.sort(ascending=False)  # in-place

    print('Saving Series as CSV')
    series.to_csv(outcsv)

    print('Unstacking Series')
    df = series.unstack()

    print('Transposing DataFrame')
    df = df.transpose()

    print('Saving DataFrame')
    df.to_pickle(outdf)


if __name__ == '__main__':
    import sys
    sort(*sys.argv[1:])
