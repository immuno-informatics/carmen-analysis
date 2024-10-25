"""Create 2D embeddings for the clustered peptides' PSSMs and group them."""

import pickle
from pathlib import Path

import pandas as pd
from sklearn.cluster import AgglomerativeClustering

project_dir = Path("/path/to/carmen-analysis")

data_dir = project_dir.joinpath("data")
data_pub_dir = data_dir.joinpath("to-be-published")
data_subs_dir = data_dir.joinpath("subsidiary-files")

samples_pssms_file = data_subs_dir.joinpath("samples-main-pssms.tsv")

umap_reducer_main_file = data_pub_dir.joinpath("umap-reducer-main.pickle")

db_main_file = data_dir.joinpath("carmen-main.parquet")

# Output files
main_samples_file = data_subs_dir.joinpath("main-output-table-1.tsv")
main_samples_file_2 = data_subs_dir.joinpath("main-output-table-2.tsv")

att_df = pd.read_csv(samples_pssms_file, sep="\t", header=0)

peptides = att_df["Peptides"]
samples = att_df["Sample"]
pssms = att_df.drop(["Peptides", "Sample"], axis=1)

# Load the UMAP model and create embeddings

with open(umap_reducer_main_file, "rb") as handle:
    reducer = pickle.load(handle)

embedding = reducer.transform(pssms)

embedding = pd.DataFrame(embedding)
embedding.columns = ["x", "y"]

# Agglomerative clustering model

cluster_number = 15  # assuming 15 clusters

cluster_model = AgglomerativeClustering(
    n_clusters=cluster_number, affinity="cosine", linkage="complete"
).fit(embedding[["x", "y"]].to_numpy())

embedding["label"] = cluster_model.labels_

# Create the two output tables

# Additionally, extract haplotype info from the main database
db_main = pd.read_parquet(
    db_main_file, columns=["Sample_name", "Haplotype"]
).drop_duplicates()
haplotypes = (
    pd.DataFrame(samples)
    .merge(db_main, how="inner", left_on="Sample", right_on="Sample_name")["Haplotype"]
    .reset_index(drop=True)
)
haplotypes.name = "Binding_alleles"

main_samples = pd.concat([samples, haplotypes, embedding], axis=1)
main_samples.to_csv(main_samples_file, sep="\t", index=False)

main_samples_2 = pd.concat([samples, peptides], axis=1)
main_samples_2.to_csv(main_samples_file_2, sep="\t", index=False)
