"""Create and save a UMAP model for GibbsCluster-driven samples' PSSMs 2D embeddings."""

import pickle
from pathlib import Path

import pandas as pd
import umap

project_dir = Path("/path/to/carmen-analysis")

data_dir = project_dir.joinpath("data")
data_pub_dir = data_dir.joinpath("to-be-published")
data_subs_dir = data_dir.joinpath("subsidiary-files")

samples_pssms_file = data_subs_dir.joinpath("samples-main-pssms.tsv")

umap_reducer_main_file = data_pub_dir.joinpath("umap-reducer-main.pickle")

att_df = pd.read_csv(samples_pssms_file, sep="\t", header=0)

pssms = att_df.drop(["Peptides", "Sample"], axis=1)

reducer = umap.UMAP(
    n_neighbors=50,
    min_dist=0.5,
    n_components=2,
    metric="cosine",
    random_state=42,
).fit(pssms.to_numpy())

with open(umap_reducer_main_file, "wb") as handle:
    pickle.dump(reducer, handle)
