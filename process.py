#!/usr/bin/env python3
"""Routines for processing data from its initial form into its useful form."""

import csv
import pickle

import pandas as pd
import numpy as np


def write_set(iterable, filename):
    with open(filename, 'w') as f:
        print('\n'.join(iterable), file=f)


def mut_patients(filename='data/tcga/file_manifest.txt'):
    """Write a list of patients from the manifest in the archive.

    The reason for doing this is that a patient might have no somatic
    mutations, and so it's useful to know whether the manifest says there are
    any more patients than are included in the MAF file.

    """
    with open(filename, 'r') as f:
        reader = csv.reader(f, dialect='excel-tab')
        index = next(reader).index('Barcode')
        samples = next(reader)[index].split('/')
        return set(s[:12] for s in samples)


def mut_patients_2(filename='data/tcga/curated-dna-sequencing.maf'):
    """Return a set of patients using the actual mutations maf file."""
    with open(filename, 'r') as f:
        reader = csv.reader(f, dialect='excel-tab')
        header = next(reader)
        sample_index = header.index('Tumor_Sample_Barcode')
        return set(row[sample_index][:12] for row in reader)


def mut_genes(filename='data/tcga/curated-dna-sequencing.maf'):
    """Return a set of genes that are mutated."""
    with open(filename, 'r') as f:
        reader = csv.reader(f, dialect='excel-tab')
        header = next(reader)
        gene_index = header.index('Hugo_Symbol')
        return set(row[gene_index] for row in reader)


def mutations(filename='data/tcga/curated-dna-sequencing.maf'):
    """Return a set of (patient, gene) pairs from the MAF file."""
    with open(filename, 'r') as f:
        reader = csv.reader(f, dialect='excel-tab')
        header = next(reader)
        gene_index = header.index('Hugo_Symbol')
        sample_index = header.index('Tumor_Sample_Barcode')
        pairs = []
        for row in reader:
            pairs.append((row[sample_index][:12], row[gene_index]))
        return set(pairs)


def write_mutations(pairset, filename):
    """Write a csv of patient,gene pairs."""
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(('Patient', 'Gene'))
        for pair in pairset:
            writer.writerow(pair)


def matrix(pairset, patients, genes):
    """Create the DataFrame from the mutations data file."""
    matrix = pd.DataFrame(np.zeros((len(patients), len(genes))),
                          index=patients, columns=genes, dtype=bool)
    for patient, gene in pairset:
        matrix.ix[patient, gene] = True
    return matrix


def expression(filename='data/cosmic/CosmicCompleteGeneExpression.tsv',
               patients=None, genes=None):
    """Return a list of COSMIC expression calls for the given patients & genes.

    The COSMIC file is huge, containing data for all TCGA patients.  We just
    need the ones that apply to our somatic mutation patients and genes, so
    this function reads that file, and returns a list of only those entries.
    However, this list could have duplicate expression values for a given
    (patient, gene) pair.

    """
    patients = patients if patients else mut_patients()
    genes = genes if genes else mut_genes()
    with open(filename) as f:
        reader = csv.reader(f, dialect='excel-tab')
        header = next(reader)
        sidx = header.index('SAMPLE_NAME')
        gidx = header.index('GENE_NAME')
        ridx = header.index('REGULATION')

        ex_vals = []
        for row in reader:
            patient = row[sidx][:12]
            gene = row[gidx]
            if patient in patients and gene in genes:
                ex_vals.append((patient, gene, row[ridx]))

    return ex_vals


def expval(string):
    """Convert a string representation of expression into a numeric value."""
    if string == 'under':
        return 0
    elif string == 'over':
        return 2
    elif string == 'normal':
        return 1
    else:
        print('ERROR: Bad expression value!')
        return 1


def deduplicate_expression(ex_vals, patientset, geneset):
    """Transform an expression list into a unique patient x gene matrix.

    Take the list of (patient, gene, expression) tuples and turn it into a
    patient x gene matrix.  Since COSMIC has some duplicate expression values
    for some patients, we simply record which patients are duplicated, and
    remove them from the matrix.

    """
    matrix = pd.DataFrame(index=patientset, columns=geneset, dtype=np.integer)
    done = 0
    created = 0
    appended = 0
    dup_pats = set()
    for patient, gene, expression in ex_vals:
        done += 1
        if not np.isnan(matrix.at[patient, gene]):
            dup_pats.add(patient)
            appended += 1
        else:
            matrix.at[patient, gene] = expval(expression)
            created += 1
        if done % 500000 == 0:
            print('Done combining %d/%d.  Created %d, dulicates %d.' %
                  (done, len(ex_vals), created, appended))

    print('Finished combining expression values into lists.')
    patientset = patientset.difference(dup_pats)
    matrix = matrix.reindex(index=patientset)
    return matrix, patientset


def restrict_mutations(pairset, patientset, geneset):
    """Restrict a list of mutation pairs.

    Only include pairs where the patient and gene are both contained in the
    given sets.

    """
    newset = set()
    for patient, gene in pairset:
        if patient in patientset and gene in geneset:
            newset.add((patient, gene))
    return newset


def main():
    """Perform all the data processing steps.

    Read somatic mutations and get the initial patient/gene set.  Then, read
    and filter the COSMIC data.  Then, restrict patients/genes to those
    included in COSMIC.  Then, remove duplicate patients in COSMIC.  Then,
    convert COSMIC data to a matrix, and save it.  Restrict somatic mutations
    to the final patient anfd gene sets, and save it as a matrix and a list.

    """
    print('Getting patient and gene pairs from somatic mutation data...')
    mut_pats = mut_patients()
    mut_gens = mut_genes()
    pairset = mutations()

    print('Reading & filtering COSMIC expression data...')
    exp_list = expression(patients=mut_pats, genes=mut_gens)

    print('Reshaping COSMIC expression data...')
    exp_pats, exp_genes, exp_values = zip(*exp_list)

    print('Filtering patient and gene set by COSMIC data.')
    patients = mut_pats.intersection(exp_pats)
    genes = mut_gens.intersection(exp_genes)

    print('Deduplicating COSMIC data...')
    expmtrx, patients = deduplicate_expression(exp_list, patients, genes)

    print('Saving final COSMIC data...')
    with open('data/expression.pickle', 'wb') as f:
        pickle.dump(expmtrx, f)

    print('Writing final patient and gene list.')
    write_set(patients, 'data/patients.txt')
    write_set(genes, 'data/genes.txt')

    print('Filtering somatic mutations by final patient and gene list...')
    newmutations = restrict_mutations(pairset, patients, genes)
    mutmtrx = matrix(newmutations, patients, genes)

    print('Writing final somatic mutation matrix and list.')
    write_mutations(newmutations, 'data/mutations.csv')
    with open('data/mutations.pickle', 'wb') as f:
        pickle.dump(mutmtrx, f)

    print('Tada!  Data is ready to process.')


if __name__ == '__main__':
    main()
