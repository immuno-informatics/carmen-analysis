# CARMEN Database Analysis

CARMEN database analysis with data and figure generation for the upcoming paper

## Setup

### Python Environment

We suggest to install [Conda](https://docs.anaconda.com/miniconda) to handle Python operations.

First, create a proper Python environment before using any of the scripts/notebooks (you can change the name `carmen-analysis` to anything else):

```bash
conda env create --name carmen-analysis --file environment.yml --yes
conda activate carmen-analysis
```

### Additional Databases

#### CARMEN Core Files

The analysis pipeline runs on certain core CARMEN database files available at [HERE](https://google.com). Please download the following files:

1. neoantigen-db-main.tsv
2. x

And subsequently move them to the *data* directory.

#### NMDP Registry Haplotype Frequencies

This pipeline is using a set of tables indicating frequencies of HLA class I alleles and haplotypes for certain populations that are a part of the [NMDP Registry Haplotype Frequencies database](https://frequency.nmdp.org).

The files used within the pipeline are as follows (names as downloaded from the website):

1. A.xlsx
2. B.xlsx
3. C.xlsx
4. A\~C\~B.xlsx

To be able to run the analysis and further process all necessary data, please download the database and put the above files into the *data/nmdp-hla-frequencies* directory.

#### MHC Motif Atlas

Partial data from the [MHC Motif Atlas class I alleles](http://mhcmotifatlas.org/class1) database is used during the analysis. Please use the "Download Data" button to access the "All Ligands" option and save the resulting file. The file *data_classI_all_peptides.txt* (name as downloaded from the website) should be moved to the *data* directory.

## How to Run

Step-by-step instruction how to go from the source CARMEN files to all the results.

## Published Data Files

The result files to be published can be found in the *data/to-be-published* directory.

### Peptide Clusters Characteristics Table

The *data/to-be-published/carmen-peptide-clusters-characteristics.parquet* file is a comma-separated table representing peptide clusters created from dividing CARMEN samples through GibbsCluster with all their relevant characteristics (binding alleles, UMAP representation, frequency of occurrence in a population, etc.).

The table columns are as follows:

- **ID**

  Peptide cluster identification code in the form of "P\_*\<X\>*", where *X* is a consecutive cluster number ("1", "2", "3", etc.).

- **Sample_name**

  Source sample identifier, as in the main CARMEN database file.

- **Peptides**

  Semicolon-separated list of peptides.

- **Binding_alleles**

  Semicolon-separated list of binding alleles associated with the peptides.

- **UMAP\_x**

  The x coordinate from the UMAP model created based on a PSSM of the clustered peptides.

- **UMAP\_y**

  The y coordinate from the UMAP model created based on a PSSM of the clustered peptides.

- **Label**

  Cluster label given by the clustering model based on UMAP coordinates ("1", "2", "3", etc.).

- **Freq\_*\<X\>*\_*\<Y\>***

  A collection of columns indicating frequency of occurrence of given MHC class I binding alleles (A or B or C) within certain populations, where *X* is a population code ("AAFA", "EURCAU", "NCHI", etc.) and *Y* is one of the three MHC class I genes indicators ("A", "B", or "C").

- **Freq_*\<X\>*\_ABC\_*\<Y\>***

  A collection of columns indicating frequency of occurrence of given triplets of MHC class I alleles (A+B+C) within certain populations, where *X* is a population code ("AAFA", "EURCAU", "NCHI", etc.) and *Y* indicates whether a specific triplet of alleles is accounted for if any of the alleles match ("any") or if all alleles match ("all").

### Peptide Clusters PSSMs

The *data/to-be-published/carmen-peptide-clusters-pssms.json* file is a supplement to the [peptide clusters characteristics table](#peptide-clusters-characteristics-table). Its structure is a hierarchical JSON representation of PSSMs for each of the clustered peptides (calculated based only on 9-mer peptides contained within). At the first level, P\_*\<X\>* represents individual peptide cluster IDs. Inside each cluster, the letters (A, C, ..., Y) represent the 20 amino acids. For each amino acid, there is a set of key-value pairs where the keys (1, 2, ..., 9) represent specific positions in the peptide sequence, and the values represent frequency for that amino acid at the corresponding position. An example of that structure is presented below.

```json
{
  "P_1":{
    "A":{
      "1":0.5,
      "2":0.011,
      "3":0.09,
      ...,
      "9":0.188
    },
    "C":{
      "1":0.06,
      "2":0.2,
      "3":0.0009,
      ...,
      "9":0.303
    },
    ...,
    "Y":{
      "1":0.5,
      "2":0.011,
      "3":0.09,
      ...,
      "9":0.188
    }
  },
  "P_2":{
    "A":{
      "1":0.23,
      ...,
      "9":0.0988
    },
    ...,
    "Y":{
      "1":0.641,
      ...,
      "9":0.1
    }
  },
  ...,
  "P_N":{
    ...
  }
}
```

An example of loading the PSSMs in Python would be as follows:

```python
import json

import pandas as pd

with open("data/to-be-published/carmen-peptide-clusters-pssms.json", "r") as handle:
    pssms = json.load(handle)

peptide_cluster_id = "P_1"
selected_pssm = pd.DataFrame(pssms[peptide_cluster_id])
```
