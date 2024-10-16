"""Set up all working directories before running the pipeline."""

from pathlib import Path

project_dir = Path("/data/teamgdansk/apalkowski/carmen-analysis")

data_dir = project_dir.joinpath("data")
data_pub_dir = data_dir.joinpath("to-be-published")
data_subs_dir = data_dir.joinpath("subsidiary-files")

figures_dir = project_dir.joinpath("figures")
figures_gen_dir = figures_dir.joinpath("script-generated")
figures_main_dir = figures_gen_dir.joinpath("main-panels")
figures_supp_dir = figures_gen_dir.joinpath("supplementary")

temp_dir = project_dir.joinpath("temp")

data_pub_dir.mkdir(parents=True, exist_ok=True)
data_subs_dir.mkdir(parents=True, exist_ok=True)
figures_main_dir.mkdir(parents=True, exist_ok=True)
figures_supp_dir.mkdir(parents=True, exist_ok=True)
temp_dir.mkdir(parents=True, exist_ok=True)
