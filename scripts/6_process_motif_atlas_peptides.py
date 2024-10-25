"""Process peptides from the MHC Motif Atlas database."""

import pickle
from pathlib import Path

# project_dir = Path("/path/to/carmen-analysis")
project_dir = Path("/data/teamgdansk/apalkowski/carmen-analysis")

data_dir = project_dir.joinpath("data")
data_pub_dir = data_dir.joinpath("to-be-published")
data_subs_dir = data_dir.joinpath("subsidiary-files")

# MHC Motif Atlas peptides, available from the project website (http://mhcmotifatlas.org)
motif_atlas_peps_src_file = data_dir.joinpath("data_classI_all_peptides.txt")
motif_atlas_peps_col_allele = "Allele"  # column name
motif_atlas_peps_col_peptide = "Peptide"  # column name

motif_atlas_peps_file = data_subs_dir.joinpath("motif-atlas-peptides.csv")

umap_reducer_main_file = data_pub_dir.joinpath("umap-reducer-main.pickle")

main_motif_atlas_emb_file = data_subs_dir.joinpath("main-embeddings_MHCMotifA_pssms.tsv")
main_motif_atlas_emb_col_allele = "Allele"  # column name
main_motif_atlas_emb_col_x = "0"  # column name
main_motif_atlas_emb_col_y = "1"  # column name
main_motif_atlas_emb_col_peps = "Peptides"  # column name
