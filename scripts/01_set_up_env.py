"""Set up all main working directories before running the pipeline."""

from pathlib import Path

# project_dir = Path("/path/to/carmen-analysis")
project_dir = Path("/data/teamgdansk/apalkowski/carmen-analysis")

data_dir = project_dir.joinpath("data")
data_pub_dir = data_dir.joinpath("to-be-published")
data_subs_dir = data_dir.joinpath("subsidiary-files")

gibbs_dir = data_subs_dir.joinpath("gibbscluster-results")

peps_per_sample_dir = data_subs_dir.joinpath("peptides-per-sample")

random_pep_dir = data_subs_dir.joinpath("random-peptides")

ml_models_dir = data_subs_dir.joinpath("ml-models")

figures_dir = project_dir.joinpath("figures")
figures_gen_dir = figures_dir.joinpath("script-generated")
figures_main_dir = figures_gen_dir.joinpath("main-panels")
figures_supp_dir = figures_gen_dir.joinpath("supplementary")

temp_dir = project_dir.joinpath("temp")

extra_dirs = [
    data_pub_dir,
    data_subs_dir,
    gibbs_dir,
    peps_per_sample_dir,
    figures_main_dir,
    figures_supp_dir,
    temp_dir,
    random_pep_dir,
    ml_models_dir,
]

for ed in extra_dirs:
    ed.mkdir(parents=True, exist_ok=True)
