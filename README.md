eecs459
=======
EECS 459 Project Code


Setup
-----

To run this, you need all of the following:

* Python 3
    * Should be installed by default on Ubuntu 14.  Run with `python3`.
* NumPy
    * `sudo apt-get install python3-numpy`
* Pandas
    * `sudo apt-get install python3-pandas`

Also you'll need the TCGA and COSMIC data mentioned below, placed in the correct
locations.

To get the data from its initial downloaded form into the form we'd like it for
analysis, I wrote the file `process.py`.  You should be able to run it by simply
issuing the command `./process.py`.


Data
----

Data is TCGA BRCA Somatic Mutation (curated), and COSMIC gene expression.

### Initial Files

Data from TCGA went into the `data/tcga` folder, and data from COSMIC went into
`data/cosmic`.  Inside the COSMIC folder, you have just the single file,
`CosmicCompleteGeneExpression.tsv`, which is gargantuan.  Inside the TCGA
folder, we have:

* `file_manifest.txt`
* `curated-dna-sequencing.maf` (this is the actual data file, renamed)

My data pre-processing program, `process.py`, expects these filenames and
locations.

### Statistics about the data

#### Somatic mutations

* Number of patients with mutation data: 977
* Number of genes with mutation data: 17689
* Number of (patient, gene) mutation pairs: 80,006
* Total number of possible pairs: 977 x 17689 = 17,282,153

#### Gene expression

I filtered the COSMIC gene expression to only include the ones with the patient
and gene we have in the somatic mutation data.

* Total (patient, gene) expression pairs: 15,000,220
* Total unique (patient, gene) expression pairs: 14,905,780
* This boils down to what appears to be 947 patients and 15,740 genes.
* The difference between the two numbers is that there are 6 patients with a
  second set of gene expression values (953 patients times 15,740 genes is
  15,000,220 pairs).
* Since these 6 patients have frequent disagreements in gene expression, we will
  exclude them, leaving 941 patients.

### Final processed data

After removing duplicates from the gene expression data, it appears we should
have the following:

* Total patients: 941
* Total genes: 15,740
* Gene expression values for all 14,811,340 pairs of patients and genes.
* Total recorded mutations: 71,695

Both the somatic mutations and the gene expression are expressed in a matrix,
with rows as patients and columns as genes.  The somatic mutations are a matrix
of booleans, while the gene expressions are a matrix of integers (0, 1, and 2
for under, normal, and over -expressed respectively).

The following files should be produced by my data processing program:

* `data/mutations.pickle` containing the somatic mutation matrix
* `data/mutations.csv` containing a list of patient,gene mutation pairs
* `data/expression.pickle` containing the gene expression matrix
* `data/patients.txt` containing a list of patients
* `data/genes.txt` containing a list of genes

### Reproducing

If you have the data in the correct places above, you should be able to simply
run `./process.py`, and you will get the expression and mutations dataframes.


Pairwise Mutual Information
---------------------------

### Math Functions

The definitions for the major functions in computing pairwise mutual information
are all in `mathfunc.py`.  I copied the implementations from my `tcga` module
and modified them as necessary.

### Experiment

I have an experiment class, which I use for doing highly parallel tasks in which
you only need to vary parameters.  I copied that class into `experiment.py`
(again, from my TCGA module).

### Mutual Information Computation

To compute pairwise mutual information, I subclassed `Experiment` to
`MiExperiment` in `mi_computation.py`.  This uses the math functions to do all
pairs mutual information.  Then, it stores the mutual information values it
computed, both into a CSV and into a DataFrame.

There's a bit of a hack in order to make this all work.  The way I stored the
results is a bit funky -- theres a `storage.py` module that has a class which
will store each batch of data as it comes in.  This is because if I didn't I
would have to copy all these DataFrames into each process, and it would have
been awful (200GB of memory, anyone?).

### Results

The results I got were ultimately:

* `data/ee.csv`: expression-expression pairs
* `data/em.csv`: expression-mutation pairs
* `data/me.csv`: mutation-expression pairs
* `data/mm.csv`: mutation-mutation pars

The process crashed before I could pickle the DataFrames.  Which actually wasn't
the worst thing to do.
