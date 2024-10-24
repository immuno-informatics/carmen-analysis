#!/bin/bash -l

project_dir="/path/to/the/main/project/dir"

input_dir=$project_dir"/data/subsidiary-files/peptides-per-sample"

sbatch_script=$project_dir"/scripts/gibbscluster-sbatch.sh"

file_list=$(find $input_dir -maxdepth 1 -type f -name "*.tsv")

for f in ${file_list[@]}; do
  f_name=$(basename $f)
  sbatch --job-name="gibbs_cl3-7-"$f_name $sbatch_script $f
done
