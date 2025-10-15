import polars
import sys
import math
import numpy as np
from pathlib import Path


project_dir = Path("/data/teamgdansk/mwaleron/carmen-analysis")
temp_dir = project_dir.joinpath("temp")
calculation_spam_dir = temp_dir.joinpath("parquetspam")
data_dir = project_dir.joinpath("data")
data_subs_dir = data_dir.joinpath("subsidiary-files")
main_samples_file_3 = data_subs_dir.joinpath("emilia_umap_with_ids.parquet")

braun_full = polars.read_parquet(data_subs_dir.joinpath("braun_mutations_with_idx.parquet"))
braun_to_score = polars.read_parquet(data_subs_dir.joinpath("narrow_and_broad_braun_mutation_wo_prom.parquet"))
emilia_umap_with_ids = polars.read_parquet(main_samples_file_3)
xallcen = (emilia_umap_with_ids["x"].min() + emilia_umap_with_ids["x"].max() )/2
yallcen = (emilia_umap_with_ids["y"].min() + emilia_umap_with_ids["y"].max() )/2

def promiscuity_peptide_list(input):
    scorebase = emilia_umap_with_ids\
        .with_columns(
        scafres = polars.col("Peptides").list.set_intersection(input
        ))\
        .with_columns(alen = polars.col("scafres").list.len()).sort("alen")
    xy = scorebase.filter(polars.col("alen")>0).select("x","y")
    count = len(xy)
    if count == 0 :
        return 0
    xcen = xy["x"].sum()/count
    ycen = xy["y"].sum()/count
    lth = np.sqrt( np.power(xallcen-xcen,2) + np.power(yallcen-ycen,2))
    popcov_but_sqrt = count / math.sqrt(1+lth)
    return popcov_but_sqrt

braun_narrow_prom = []
braun_broad_prom = []
for row in braun_to_score.iter_rows():
    braun_narrow_prom.append(promiscuity_peptide_list(row[1]))
    braun_broad_prom.append(promiscuity_peptide_list(row[10]))  
    
res=braun_to_score.hstack(polars.DataFrame({"Promiscuity_narrow":braun_narrow_prom, "Promiscuity_broad":braun_broad_prom}))
res.write_parquet("braun_prom_just_overlap.parquet")
complete_res = res.join(braun_full, on="mutid")
complete_res.write_parquet("braun_prom_all.parquet")