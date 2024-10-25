# CARMEN Database Analysis

CARMEN database analysis with data and figure generation for the [upcoming](https://www.google.com) paper

## Setup

**ADD INFO ABOUT CLONING ETC.**

Some of the pipeline files and file processing descriptions included assume work in a Linux environment.

### Python Environment

We suggest to install [Conda](https://docs.anaconda.com/miniconda) to handle Python operations.

The described pipeline used multiple Python environments to handle different parts of the analysis while preparing the results. Details of those environments have been saved as [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) environment specifications to `environment-\*.yml` files. Please follow in detail the instructions included in the [How to Run the Analysis](#how-to-run-the-analysis) section to achieve reproducibility.

### Required Databases

The analysis pipeline requires certain external files to produce all results. Please follow the subsections below for details.

#### CARMEN Core Files

The analysis pipeline runs on certain core CARMEN database files available at [the main repository](https://doi.org/10.5281/zenodo.13928442). Please download the following files:

1. `carmen-main.parquet`&mdash;the main part of the database. Contains a dataset gathered from standardized reprocessing of 72 publicly available immunopeptidomic mass spectrometry datasets. Contains all gathered peptides and their annotations.
2. x

And subsequently move them to the *data* directory.

#### CARMEN Analysis Additional Files

Certain aspects of the analysis and its products posses a degree of stochasticity that could alter the results in a noticeable way when run again. Because of that a number of files which are necessary for subsequent steps of the analysis were made available at the [CARMEN analysis paper repository](https://google.com) for reproducibility. You may download the following files to reproduce the analysis:

1. `umap-reducer-main.pickle`&mdash;serialized Python object that represents a trained UMAP from the `umap-learn` package.
2. x

And subsequently move them to the `data/to-be-published` directory.

#### NMDP Registry Haplotype Frequencies

This pipeline is using a set of tables indicating frequencies of HLA class I alleles and haplotypes for certain populations that are a part of the [NMDP Registry Haplotype Frequencies database](https://frequency.nmdp.org).

The files used within the pipeline are as follows (names as downloaded from the website):

1. `A.xlsx`
2. `B.xlsx`
3. `C.xlsx`
4. `A\~C\~B.xlsx`

To be able to run the analysis and further process all necessary data, please download the database and put the above files into the `data/nmdp-hla-frequencies` directory.

#### MHC Motif Atlas

Partial data from the [MHC Motif Atlas class I alleles](http://mhcmotifatlas.org/class1) database is used during the analysis. Please use the "Download Data" button to access the "All Ligands" option and save the resulting file. The file `data_classI_all_peptides.txt` (name as downloaded from the website) should be moved to the `data` directory.

## How to Run the Analysis

To be able to run the entire analysis and produce its results (data files, figures, etc.) we need to set-up and run specific scripts in a particular order. Please follow the instructions below.

**Always remember to review all paths in the scripts and make changes appropriate to your own system's environment.**

### 1. Create Working Directories

First, create a new Python environment based on a proper specification (`environment-1-6.yml`) before using any of the scripts/notebooks described below (you can change the name `carmen-analysis-1-6` to anything else):

```bash
conda env create --name carmen-analysis-1-6 --file environment-1-6.yml --yes
conda activate carmen-analysis-1-6
```

After which, run the `scripts/1_set_up_env.py` script:

```bash
python scripts/1_set_up_env.py
```

### 2. Prepare Samples for GibbsCluster

Run `scripts/2_prep_gibbscluster_samples.py` script:

```bash
python scripts/2_prep_gibbscluster_samples.py
```

### 3. Run GibbsCluster Over the Samples

Cluster peptides within samples using [GibbsCluster-2.0](https://services.healthtech.dtu.dk/services/GibbsCluster-2.0).

Because of the big number of samples and peptides to analyze, here we utilized a computing cluster with the [Slurm Workload Manager](https://slurm.schedmd.com) system.

There are two script files used to produce [GibbsCluster](https://services.healthtech.dtu.dk/services/GibbsCluster-2.0) results:

- `scripts/gibbscluster-sbatch.sh`&mdash;sbatch script with configuration for a singular Slurm job of running [GibbsCluster](https://services.healthtech.dtu.dk/services/GibbsCluster-2.0) with given arguments
- `scripts/3_gibbscluster-spawn-jobs.sh`&mdash;Bash script that iterates over input files and, for each, creates a Slurm job to run it through [GibbsCluster](https://services.healthtech.dtu.dk/services/GibbsCluster-2.0)

The entry point is the `scripts/3_gibbscluster-spawn-jobs.sh` script where you set up all necessary input and output configuration. Remember to set up a proper job configuration in `scripts/gibbscluster-sbatch.sh`, appropriate for your [Slurm](https://slurm.schedmd.com) system configuration, and then run the main script:

```bash
bash scripts/3_gibbscluster-spawn-jobs.sh
```

In case there is no access to a [Slurm](https://slurm.schedmd.com) system, you can try to run [GibbsCluster](https://services.healthtech.dtu.dk/services/GibbsCluster-2.0) serially (or adapt a similar way):

```bash
project_dir="/path/to/carmen-analysis"

input_dir=$project_dir"/data/subsidiary-files/peptides-per-sample"

output_dir=$project_dir"/data/subsidiary-files/gibbscluster-results"

file_list=$(find $input_dir -maxdepth 1 -type f -name "*.tsv")

for f in ${file_list[@]}; do
  f_name=$(basename $f)
  output_file=$output_dir"/output_file_"$f_name".csv.out"
  gibbscluster -f $f -R $output_dir -P $f_name".out-len-9" -T -j 10 -g3-7 -l9 -S1 > $output_file
done
```

### 4. Extract PSSMs of Clustered Peptides

Extract position-specific scoring matrices from [GibbsCluster](https://services.healthtech.dtu.dk/services/GibbsCluster-2.0) results. Run `scripts/4_extract_pssms.py` script:

```bash
python scripts/4_extract_pssms.py
```

### 5. Create Embeddings for PSSMs and Group Them

Take the position-specific scoring matrices calculated for clustered peptides from samples and create 2-dimensional embeddings using a UMAP model. Next, group those representations using a clustering algorithm.

To use the UMAP model trained and saved by the authors, follow the instructions in the [CARMEN Analysis Additional Files](#carmen-analysis-additional-files) section for the `umap-reducer-main.pickle` file, and then run the `scripts/5_umap_embeddings_and_labels.py` script:

```bash
python scripts/5_umap_embeddings_and_labels.py
```

The exact way how to re-create (or re-train) the UMAP model is contained in the `scripts/umap_model_main.py` script. Running this script will overwrite any saved model. Please be advised that recreating the UMAP model on a different machine setup may alter a significant portion of the results:

```bash
python scripts/umap_model_main.py
```

After creating a new model run the `scripts/5_umap_embeddings_and_labels.py` script to calculate new embeddings and labels.

### 6. Process the MHC Motif Atlas Peptides

Filter-out non-human alleles and peptides of lengths other than 9 amino acids, and then create 2-dimensional embeddings for each unique allele by using the UMAP model that was previously trained on peptides from CARMEN samples. Run the `scripts/6_process_motif_atlas_peptides.py` script:

```bash
python scripts/6_process_motif_atlas_peptides.py
```

### 7. OMG

First, create a new Python environment based on a proper specification (`environment-7-X.yml`) before using any of the scripts/notebooks described below (you can change the name `carmen-analysis-7-X` to anything else); deactivate the previous environment if you have it loaded:

```bash
# Run this line if you have another Conda environment loaded:
conda deactivate

conda env create --name carmen-analysis-7-X --file environment-7-X.yml --yes
conda activate carmen-analysis-7-X
```

### X. Asd

Qwe.

## Published Data Files

The result files to be published can be found in the `data/to-be-published` directory.

### Peptide Clusters Characteristics Table

The `data/to-be-published/carmen-peptide-clusters-characteristics.parquet` file is an [Apache Parquet](https://parquet.apache.org) data file table representing peptide clusters created from dividing CARMEN samples through [GibbsCluster](https://services.healthtech.dtu.dk/services/GibbsCluster-2.0) with all their relevant characteristics (binding alleles, UMAP representation, frequency of occurrence in a population, etc.). Key column: **ID**.

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

The `data/to-be-published/carmen-peptide-clusters-pssms.json` file is a supplement to the [peptide clusters characteristics table](#peptide-clusters-characteristics-table). Its structure is a hierarchical JSON representation of PSSMs for each of the clustered peptides (calculated based only on 9-mer peptides contained within). At the first level, P\_*\<X\>* represents individual peptide cluster IDs. Inside each cluster, the letters (A, C, ..., Y) represent the 20 amino acids. For each amino acid, there is a set of key-value pairs where the keys (1, 2, ..., 9) represent specific positions in the peptide sequence, and the values represent frequency for that amino acid at the corresponding position. An example of that structure is presented below.

```json
{
  "P_1":{
    "A":{
      "1":0.5,
      "2":0.011,
      "3":0.09,
      //...,
      "9":0.188
    },
    "C":{
      "1":0.06,
      "2":0.2,
      "3":0.0009,
      //...,
      "9":0.303
    },
    //...,
    "Y":{
      "1":0.5,
      "2":0.011,
      "3":0.09,
      //...,
      "9":0.188
    }
  },
  "P_2":{
    "A":{
      "1":0.23,
      //...,
      "9":0.0988
    },
    //...,
    "Y":{
      "1":0.641,
      //...,
      "9":0.1
    }
  },
  //...,
  "P_X":{
    //...
  }
}
```

An example of loading the PSSMs in Python would be as follows:

```python
import json

import pandas as pd

with open("data/to-be-published/carmen-peptide-clusters-pssms.json", "r") as handle:
    pssms = json.load(handle)

all_peptide_cluster_ids = pssms.keys()
peptide_cluster_id = "P_1"
selected_pssm = pd.DataFrame(pssms[peptide_cluster_id])
```
