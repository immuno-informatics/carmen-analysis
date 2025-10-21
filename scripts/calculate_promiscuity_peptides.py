import polars as pl
import numpy as np
import math
import sys
from pathlib import Path

project_dir = Path("/data/teamgdansk/mwaleron/carmen-analysis")
data_dir = project_dir.joinpath("data")
data_subs_dir = data_dir.joinpath("subsidiary-files")
main_samples_file_3 = data_subs_dir.joinpath("emilia_umap_with_ids.parquet")
emilia_umap_with_ids = pl.read_parquet(main_samples_file_3)
temp_dir = project_dir.joinpath("temp")


i = sys.argv[1]
xallcen = (emilia_umap_with_ids["x"].min() + emilia_umap_with_ids["x"].max() )/2
yallcen = (emilia_umap_with_ids["y"].min() + emilia_umap_with_ids["y"].max() )/2

def promiscuity_single_peptide(input):
    scorebase = emilia_umap_with_ids.with_columns(hit = pl.col("Peptides").list.contains(input))
    xy = scorebase.filter(pl.col("hit")==True).select("x","y")
    count = len(xy)
    if count == 0 :
        return 0
    xcen = xy["x"].sum()/count
    ycen = xy["y"].sum()/count
    lth = np.sqrt( np.power(xallcen-xcen,2) + np.power(yallcen-ycen,2))
    popcov_but_sqrt = count / math.sqrt(1+lth)
    return popcov_but_sqrt

znow_bijatyka = pl.read_parquet(temp_dir.joinpath(f"trashpanda{i}.parquet"))
znow_bijatyka.with_columns(
    score = pl.col("Peptide").map_elements(
        promiscuity_single_peptide, 
        return_dtype=pl.Float64()))\
    .write_parquet(temp_dir.joinpath(f"bijatyka_caly_dzien{i}.parquet"))
