"""Storage of data."""

import numpy as np
import pandas as pd

_eecsv = 'data/ee.csv'
_emcsv = 'data/em.csv'
_mecsv = 'data/me.csv'
_mmcsv = 'data/mm.csv'


class Storage:
    def __init__(self, genes):
        self._ee = pd.DataFrame(index=genes, columns=genes, dtype=float)
        self._em = pd.DataFrame(index=genes, columns=genes, dtype=float)
        self._me = pd.DataFrame(index=genes, columns=genes, dtype=float)
        self._mm = pd.DataFrame(index=genes, columns=genes, dtype=float)

    def append(self, fn, series):
        gene = series.name
        with open(fn, 'a') as f:
            for geneB in series.index:
                if not np.isnan(series[geneB]):
                    print('%s, %s, %f' % (gene, geneB, series[geneB]), file=f)

    def store_ee(self, series):
        self._ee[series.name] = series
        self.append(_eecsv, series)

    def store_em(self, series):
        self._em[series.name] = series
        self.append(_emcsv, series)

    def store_me(self, series):
        self._me[series.name] = series
        self.append(_mecsv, series)

    def store_mm(self, series):
        self._mm[series.name] = series
        self.append(_mmcsv, series)
