"""Mutual information computation!

Contains an Experiment that will compute the mutual information between every
pair of genes, mutation and expression.  It uses multiprocessing to efficiently
distribute the load across CPU cores.

"""

import pickle

import pandas as pd
import numpy as np

from experiment import Experiment
import mathfunc as mf
from storage import Storage

storage = None


class MiExperiment(Experiment):

    def __init__(self):
        global storage
        super().__init__()
        # Open the data files.
        print('Opening expression and mutation pickles...')
        with open('data/expression.pickle', 'rb') as f:
            self.expression = pickle.load(f)
        with open('data/mutations.pickle', 'rb') as f:
            self.mutations = pickle.load(f)

        print('Initializing the storage...')
        storage = Storage(self.expression.columns)

        print('Precomputing entropy...')
        # Precompute the entropy for the two matrices.
        self.expression_entropy = mf.precompute_entropy(self.expression,
                                                        domain=(0, 1, 2))
        self.mutation_entropy = mf.precompute_entropy(self.mutations,
                                                      domain=(0, 1))
        self._params['gene'] = list(reversed(self.expression.columns))

    def task(self, config):
        gene = config[0]
        return mf.all_pairs_mutual_info(self.expression, self.mutations,
                                        self.expression_entropy,
                                        self.mutation_entropy, gene)

    def result(self, retval):
        ee, em, me, mm = retval
        storage.store_ee(ee)
        storage.store_em(em)
        storage.store_me(me)
        storage.store_mm(mm)

    def append_series(self, gene, fname, series):
        with open(fname, 'a') as f:
            for geneB in series.index:
                if not np.isnan(series[geneB]):
                    print('%s, %s, %f' % (gene, geneB, series[geneB]), file=f)
