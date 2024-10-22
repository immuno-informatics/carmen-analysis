"""."""

from pathlib import Path

import pandas as pd

project_dir = Path("/data/teamgdansk/apalkowski/carmen-analysis")

data_dir = project_dir.joinpath("data")
data_subs_dir = data_dir.joinpath("subsidiary-files")

gibbs_dir = data_subs_dir.joinpath("gibbscluster")
peps_per_sample_dir = data_subs_dir.joinpath("peptides-per-sample")

# data_gibbs_dir.mkdir(parents=True, exist_ok=True)
# samples_peptides_dir.mkdir(parents=True, exist_ok=True)

# split the files

df = pd.read_csv("neo_uniq_sample_peptide_list.csv", sep="\t")

grouped_by_sample = df.groupby(["Sample_name"])
unique_keys = df["Sample_name"].unique()

for unique_key in unique_keys:
    grouped_by_sample.get_group(unique_key).to_csv(
        peps_per_sample_dir.joinpath(unique_key, ".csv"), index=False
    )
