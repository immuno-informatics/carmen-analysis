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
dataslice = polars.read_parquet(calculation_spam_dir.joinpath(f"scaf20_slice_{i}.parquet"))
emilia_umap_with_ids = polars.read_parquet(main_samples_file_3)
def vector_properties(x1, y1, x2, y2):
    length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    return length
xallcen = (emilia_umap_with_ids["x"].min() + emilia_umap_with_ids["x"].max() )/2
yallcen = (emilia_umap_with_ids["y"].min() + emilia_umap_with_ids["y"].max() )/2

popcov_but_sqrt = []
popcov_but_sqrt5 = []
popcov_but_sqrt10 = []
popcov_but_sqrt50 = []
popcov_but_sqrt100 = []
for scaffold in dataslice.iter_rows(named=True):
    tmp = emilia_umap_with_ids.with_columns(scafres = polars.col("Peptides").list.set_intersection(scaffold["pep_list"])).with_columns(alen = polars.col("scafres").list.len()).sort("alen")
    x = tmp.select(polars.col("x").filter(polars.col("alen")>0))["x"]
    y = tmp.select(polars.col("y").filter(polars.col("alen")>0))["y"]
    x2 = tmp.select(polars.col("x").filter(polars.col("alen")>1))["x"]
    y2 = tmp.select(polars.col("y").filter(polars.col("alen")>1))["y"]
    x3 = tmp.select(polars.col("x").filter(polars.col("alen")>2))["x"]
    y3 = tmp.select(polars.col("y").filter(polars.col("alen")>2))["y"]
    x4 = tmp.select(polars.col("x").filter(polars.col("alen")>3))["x"]
    y4 = tmp.select(polars.col("y").filter(polars.col("alen")>3))["y"]
    
    count = len(x)
    count5 = len(x2)
    count10 = len(x3)
    count50 = len(x4)
    count = len(x)
    if count == 0:
        popcov_but_sqrt.append(0)
        popcov_but_sqrt5.append(0)
        popcov_but_sqrt10.append(0)
        popcov_but_sqrt50.append(0)
        continue
    xcen = x.sum()/len(x)
    ycen = y.sum()/len(y)
    length = vector_properties(xcen, ycen, xallcen, yallcen)
    popcov_but_sqrt.append(count / math.sqrt(1+length))
    
    
    if count5 == 0:
        popcov_but_sqrt5.append(0)
        popcov_but_sqrt10.append(0)
        popcov_but_sqrt50.append(0)
        continue
    
    x2cen = x2.sum()/len(x2)
    y2cen = y2.sum()/len(y2)
    length5 = vector_properties(x2cen, y2cen, xallcen, yallcen)
    popcov_but_sqrt5.append(count5 / math.sqrt(1+length5))
    
    if count10 == 0:
        popcov_but_sqrt10.append(0)
        popcov_but_sqrt50.append(0)
        continue
    
    x10cen = x3.sum()/len(x3)
    y10cen = y3.sum()/len(y3)
    length10 = vector_properties(x10cen, y10cen, xallcen, yallcen)
    popcov_but_sqrt10.append(count10 / math.sqrt(1+length10))
    
    if count50 == 0:
        popcov_but_sqrt50.append(0)
        continue
    
    x50cen = x4.sum()/len(x4)
    y50cen = y4.sum()/len(y4)
    length50 = vector_properties(x50cen, y50cen, xallcen, yallcen)
    popcov_but_sqrt50.append(count50 / math.sqrt(1+length50))
    
res=dataslice.hstack(polars.DataFrame({"popcov_but_sqrt":popcov_but_sqrt, "popcov_but_sqrt2":popcov_but_sqrt5, "popcov_but_sqrt3":popcov_but_sqrt10, "popcov_but_sqrt4":popcov_but_sqrt50}))
res.write_parquet(calculation_spam_dir.joinpath(f"rescored_scaf20_slice_{i}.parquet"))