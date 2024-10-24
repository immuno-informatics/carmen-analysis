"""."""

from pathlib import Path

import pandas as pd

project_dir = Path("/data/teamgdansk/apalkowski/carmen-analysis")

data_dir = project_dir.joinpath("data")
data_subs_dir = data_dir.joinpath("subsidiary-files")

peps_per_sample_dir = data_subs_dir.joinpath("peptides-per-sample")

db_main_file = data_dir.joinpath("carmen-main.parquet")

# Extract "Sample_name" and "Peptide" columns from the main DB

db_main = pd.read_parquet(db_main_file)

pep_sample = db_main[["Sample_name", "Peptide"]]
pep_sample = (
    pep_sample.groupby("Peptide")["Sample_name"]
    .apply(lambda x: ",".join(x))
    .reset_index()
)

del db_main

# Save .tsv files with a peptide list per each sample

# Split the 'Sample_name' column by commas, then explode it into individual rows
pep_sample["Sample_name"] = pep_sample["Sample_name"].str.split(",")
pep_sample = pep_sample.explode("Sample_name")

# Clean any leading/trailing spaces from 'Sample_name'
pep_sample["Sample_name"] = pep_sample["Sample_name"].str.strip()

# Group the DataFrame by 'Sample_name'
grouped_by_sample = pep_sample.groupby("Sample_name")

# Iterate over each unique 'Sample_name' and save the corresponding peptides
for sample_name, group in grouped_by_sample:
    # Clean up the sample name to ensure it's a valid file name
    cleaned_sample_name = sample_name.replace("/", "_").replace(" ", "_")

    # Create a TSV file for each sample, only containing the 'Peptide' column
    output_file = peps_per_sample_dir.joinpath(f"{cleaned_sample_name}.tsv")

    # Save only the 'Peptide' column for the current sample
    group[["Peptide"]].to_csv(output_file, sep="\t", header=False, index=False)
