"""Math functions for EECS 459 Project.

Implementations for entropy and mutual information of discrete datasets.

Disclaimer: Some of this code was taken from my TCGA module:
https://bitbucket.org/brenns10/tcga (private)

"""

import numpy as np
import pandas as pd


def entropy(ds, domain=(0, 1)):
    """
    Computes the entropy (in bits) of a dataset.

    :param ds: The dataset to compute the entropy of.
    :param domain: The domain of the dataset.
    :return: The entropy of the dataset, in bits.
    """
    currentropy = 0
    total = len(ds)
    for value in domain:
        probability = len(ds[ds == value]) / total
        if probability != 0:  # avoid log(0)
            currentropy -= probability * np.log2(probability)

    return currentropy


def mutual_info(ds1, ds2, e1, e2, ds1domain=2, ds2domain=3):
    """Compute mutual information between two datasets.

    Uses precomputed entropy values for ds1 and ds2, so that you can cache them
    for better performance (they're reused a lot).  The default domains are set
    for ds1 to be mutations [0,1], and ds2 to be expression [0,2].

    """
    combined = ds1 + ds2 * ds1domain
    return e1 + e2 - entropy(combined, domain=range(ds1domain * ds2domain))


def precompute_entropy(matrix, domain):
    """Compute the entropy for every gene in the matrix.

    It turns out that you need to compute the entropy for a single gene a lot.
    I figured that you could precompute them all in advance and save a lot of
    time later.  And, it turns out that once you've done that, the time to
    compute mutual info is reduced by over 25% (according to some very
    unscientific timing results).

    """
    entropies = pd.Series(index=matrix.columns)
    for gene in matrix:
        entropies[gene] = entropy(matrix[gene], domain)
    return entropies


def pairwise_mutual_info(expression, mutations, expression_entropy,
                         mutation_entropy, geneA, geneB):
    """Compute all four mutual informations between two genes.

    Params: the two matrices, the entropy caches, and the two genes.
    Returns: expression-expression mutual info
             expression-mutation   mutual info
             mutation  -expression mutual info
             mutation  -mutation   mutual info.
    (unless the two genes are the same.  then it just computes expression-
    mutation mutual info).

    """

    am = mutations[geneA]
    bm = mutations[geneB]
    ame = mutation_entropy[geneA]
    bme = mutation_entropy[geneB]
    ae = expression[geneA]
    be = expression[geneB]
    aee = expression_entropy[geneA]
    bee = expression_entropy[geneB]

    ae_bm = mutual_info(ae, bm, aee, bme, 3, 2)
    if geneA == geneB:
        return ae_bm

    ae_be = mutual_info(ae, be, aee, bee, 3, 3)
    am_bm = mutual_info(am, bm, ame, bme, 2, 2)
    am_be = mutual_info(am, be, ame, bee, 2, 3)
    return ae_be, ae_bm, am_be, am_bm


def all_pairs_mutual_info(expression, mutations, expression_entropy,
                          mutation_entropy, gene):
    """Compute the mutual information for all pairs starting with gene.

    Start with the first column in the expression dataframe, and compute all
    the mutual information pairs until you reach the gene itself.  Then stop.

    Returns four Series, in the same order as pairwise_mutual_info().  They are
    each named after this gene, so that they can be concatenated in directly
    with a result dataframe.

    """
    ee = pd.Series(name=gene, index=expression.columns)
    em = pd.Series(name=gene, index=expression.columns)
    me = pd.Series(name=gene, index=expression.columns)
    mm = pd.Series(name=gene, index=expression.columns)

    geneAidx = np.where(expression.columns == gene)[0][0]
    done = 0
    for geneBidx in range(geneAidx):
        geneB = expression.columns[geneBidx]
        res = pairwise_mutual_info(expression, mutations, expression_entropy,
                                   mutation_entropy, gene, geneB)
        ee[geneB] = res[0]
        em[geneB] = res[1]
        me[geneB] = res[2]
        mm[geneB] = res[3]
        done += 1
        if done % 1000 == 0:
            print(done)
    final = pairwise_mutual_info(expression, mutations, expression_entropy,
                                 mutation_entropy, gene, gene)
    em[gene] = final
    me[gene] = final
    return ee, em, me, mm
