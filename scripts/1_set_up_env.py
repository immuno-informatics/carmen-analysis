"""Set up all working directories before running the pipeline."""

from pathlib import Path

project_dir = Path("/data/teamgdansk/apalkowski/carmen-analysis")

data_dir = project_dir.joinpath("data")
data_pub_dir = data_dir.joinpath("to-be-published")
data_subs_dir = data_dir.joinpath("subsidiary-files")

data_gibbs_dir = data_subs_dir.joinpath("gibbscluster")

samples_peptides_dir = data_subs_dir.joinpath("samples-peptides")

figures_dir = project_dir.joinpath("figures")
figures_gen_dir = figures_dir.joinpath("script-generated")
figures_main_dir = figures_gen_dir.joinpath("main-panels")
figures_supp_dir = figures_gen_dir.joinpath("supplementary")

temp_dir = project_dir.joinpath("temp")

extra_dirs = [
    data_pub_dir,
    data_subs_dir,
    data_gibbs_dir,
    samples_peptides_dir,
    figures_main_dir,
    figures_supp_dir,
    temp_dir,
]

for ed in extra_dirs:
    ed.mkdir(parents=True, exist_ok=True)
