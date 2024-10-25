"""
Iterate over the GibbsCluster-driven samples directories (n=2323) and create
a single .tsv file containing sample names, corresponding peptides, and PSSM matrices.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

project_dir = Path("/path/to/carmen-analysis")

data_dir = project_dir.joinpath("data")
data_subs_dir = data_dir.joinpath("subsidiary-files")

gibbs_dir = data_subs_dir.joinpath("gibbscluster-results")

samples_pssms_file = data_subs_dir.joinpath("samples-main-pssms.tsv")

samples = []
matrices = []
peptides = []
kld_values = []

# List sample directories
sample_dirs = [
    d for d in gibbs_dir.iterdir() if d.is_dir() and d.name.startswith("sample_")
]

for s_dir in tqdm(sample_dirs):
    sample_dir_name = s_dir.name
    sample_name = sample_dir_name.split(".tsv.")[0]

    cluster_scores_file = s_dir.joinpath("images", "gibbs.KLDvsClusters.tab")
    try:
        kld = pd.read_csv(cluster_scores_file, sep="\t")
    except FileNotFoundError:
        continue

    # max kld
    kld_max_index_number = kld.sum(axis=1).idxmax()

    if kld_max_index_number == 0:
        continue

    res_file = s_dir.joinpath("res", f"gibbs.{kld_max_index_number}g.ds.out")
    res = pd.read_csv(res_file, sep=r"\s+", usecols=["Gn", "Sequence"])

    kld_list = int(res_file.stem.split(".")[1][0])
    kld_values.append(kld_list)

    matrices_dir = s_dir.joinpath("matrices")

    mat_files = [
        f
        for f in matrices_dir.iterdir()
        if f.is_file() and f.name.endswith(f"of{kld_max_index_number}.mat")
    ]

    for file in mat_files:
        samples.append(sample_name)
        group_no = int(file.stem.split(".")[1][0])

        df = pd.read_csv(file, sep=" ", skiprows=2, header=None)
        matrix = df.iloc[:, 2:].to_numpy().flatten().tolist()
        matrices.append(matrix)

        seq_list = res[res["Gn"] == group_no - 1]["Sequence"].to_list()
        seq_list = sorted(list(set(seq_list)))
        seq_list_str = ",".join(seq_list)
        peptides.append(seq_list_str)

matrices = np.array(matrices)

output_table = pd.DataFrame({"Sample": samples, "Peptides": peptides})

for i in range(matrices.shape[1]):
    output_table[f"PSSM_{i+1}"] = matrices[:, i]

output_table.to_csv(samples_pssms_file, sep="\t", index=False)
