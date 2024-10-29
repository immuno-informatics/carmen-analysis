import polars
import sys
import math
from pathlib import Path


project_dir = Path("/data/teamgdansk/mwaleron/carmen-analysis")
temp_dir = project_dir.joinpath("temp")
calculation_spam_dir = temp_dir.joinpath("parquetspam")
data_dir = project_dir.joinpath("data")
data_subs_dir = data_dir.joinpath("subsidiary-files")
main_samples_file_3 = data_subs_dir.joinpath("emilia_umap_with_ids.parquet")

i = sys.argv[1]
dataslice = polars.read_parquet(calculation_spam_dir.joinpath(f"contig_slice_{i}.parquet"))
emilia_umap_with_ids = polars.read_parquet(main_samples_file_3)
def vector_properties(x1, y1, x2, y2):
    length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    return length
xallcen = (emilia_umap_with_ids["x"].min() + emilia_umap_with_ids["x"].max() )/2
yallcen = (emilia_umap_with_ids["y"].min() + emilia_umap_with_ids["y"].max() )/2

popcov = []
popcov_but_sqrt = []
vlen = []
raise NotImplementedError("")
for scaffold in dataslice.iter_rows(named=True):
    tmp = emilia_umap_with_ids.with_columns(scafres = polars.col("Peptides").list.set_intersection(scaffold["pep_list"])).with_columns(alen = polars.col("scafres").list.len()).sort("alen")
    x = tmp.select(polars.col("x").filter(polars.col("alen")>0))["x"]
    y = tmp.select(polars.col("y").filter(polars.col("alen")>0))["y"]
    count = len(x)
    if count == 0:
        popcov.append(0)
        vlen.append(0)
        popcov_but_sqrt.append(0)
        continue
    xcen = x.sum()/len(x)
    ycen = y.sum()/len(y)
    length = vector_properties(xcen, ycen, xallcen, yallcen)
    popcov.append(count/(1+length))
    popcov_but_sqrt.append(count / math.sqrt(1+length))
    vlen.append(length)
res=dataslice.hstack(polars.DataFrame({"pop_cov":popcov, "vlen":vlen, "popcov_but_sqrt":popcov_but_sqrt}))
res.write_parquet(calculation_spam_dir.joinpath(f"scored_contig_slice_{i}.parquet"))