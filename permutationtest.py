#!/usr/bin/env python3
"""Permutation tests for mutual information cutoffs."""

import numpy as np
import pandas as pd

from mathfunc import entropy, mutual_info


class DiscreteRandomVariable:
    """Allows random sampling of any discrete random variable."""

    def __init__(self, domain, weights):

        self._domain = np.array(domain)
        self._weights = np.array(weights) / np.sum(weights)
        self._cumsum = np.cumsum(self._weights)

    def sample(self, amount):
        """Get `amount` random values from the distribution."""
        s = pd.Series(np.random.sample(size=amount))
        snew = s.copy()

        for value, cumsum in zip(self._domain, self._cumsum):
            snew[s <= cumsum] = value
            s[s <= cumsum] = 2  # reassign so we don't mess it up.

        return snew


def random_mutual_information(drv1, drv2, size=1000):
    """Sample the random varibles and return their mutual information."""
    a = drv1.sample(size)
    b = drv2.sample(size)

    ea = entropy(a, domain=drv1._domain)
    eb = entropy(b, domain=drv2._domain)

    return mutual_info(a, b, ea, eb, len(drv1._domain), len(drv2._domain))


def get_distributions(expressionfname):
    """Get the average distribution of a dataframe."""
    expression = pd.read_pickle(expressionfname)
    counts = {}
    for col in expression.columns:
        vc = dict(expression[col].value_counts())
        for v, c in vc.items():
            l = counts.get(v, [])
            l.append(c)
            counts[v] = l
    meancounts = {v: np.mean(clist) for v, clist in counts.items()}
    return DiscreteRandomVariable(*zip(*meancounts.items()))


def me_midist(samplesize=1000, trials=10000):
    """Get the mutual information distribution for M-E pairs."""
    mdist = get_distributions('data/mutations.pickle')
    edist = get_distributions('data/expression.pickle')
    rmi = []
    for _ in range(trials):
        rmi.append(random_mutual_information(mdist, edist, samplesize))
    return np.array(rmi)


def uniform_midist(samplesize=1000, trials=10000):
    """Get the mutual information distribution for uniform M-E pairs."""
    mdist = DiscreteRandomVariable(domain=(0, 1), weights=(1, 1))
    edist = DiscreteRandomVariable(domain=(0, 1, 2), weights=(1, 1, 1))
    rmi = []
    for _ in range(trials):
        rmi.append(random_mutual_information(mdist, edist, samplesize))
    return np.array(rmi)
