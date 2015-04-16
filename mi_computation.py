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


class MiExperiment(Experiment):

    def __init__(self):
        super().__init__()

        # CSV files to append after every trial finishes, in case of problems.
        self.csvee = 'data/ee.csv'
        self.csvmm = 'data/mm.csv'
        self.csvem = 'data/em.csv'

        # Pickles to save
        self.eepickle = 'data/ee.pickle'
        self.empickle = 'data/em.pickle'
        self.mepickle = 'data/me.pickle'
        self.mmpickle = 'data/mm.pickle'

        # Open the data files.
        with open('data/expression.pickle', 'rb') as f:
            self.expression = pickle.load(f)
        with open('data/mutations.pickle', 'rb') as f:
            self.mutations = pickle.load(f)

        # Precompute the entropy for the two matrices.
        self.expression_entropy = mf.precompute_entropy(self.expression,
                                                        domain=(0, 1, 2))
        self.mutation_entropy = mf.precompute_entropy(self.mutations,
                                                      domain=(0, 1))

        # Dataframes for all the results.
        self.eedf = pd.DataFrame(index=self.expression.columns,
                                 columns=self.expression.columns,
                                 dtype=float)
        self.emdf = pd.DataFrame(index=self.expression.columns,
                                 columns=self.expression.columns,
                                 dtype=float)
        self.medf = pd.DataFrame(index=self.expression.columns,
                                 columns=self.expression.columns,
                                 dtype=float)
        self.mmdf = pd.DataFrame(index=self.expression.columns,
                                 columns=self.expression.columns,
                                 dtype=float)

        self._params['gene'] = list(self.expression.columns)

    def task(self, config):
        gene = config[0]
        return mf.all_pairs_mutual_info(self.expression, self.mutations,
                                        self.expression_entropy,
                                        self.mutation_entropy, gene)

    def result(self, retval):
        ee, em, me, mm = retval

        # Store in data frames.
        self.eedf[ee.name] = ee
        self.emdf[em.name] = em
        self.medf[me.name] = me
        self.mmdf[mm.name] = mm

        # Write out to CSVs so we have immediate storage as well.
        self.append_series(ee.name, self.csvee, ee)
        self.append_series(em.name, self.csvem, em)
        self.append_series(me.name, self.csvem, me)
        self.append_series(mm.name, self.csvmm, mm)

    def append_series(self, gene, fname, series):
        with open(fname, 'a') as f:
            for geneB in series.index:
                if not np.isnan(series[geneB]):
                    print('%s, %s, %f' % (gene, geneB, series[geneB]), file=f)
