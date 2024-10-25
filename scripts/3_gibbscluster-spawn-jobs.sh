#!/bin/bash -l

project_dir="/path/to/carmen-analysis"

input_dir=$project_dir"/data/subsidiary-files/peptides-per-sample"

output_dir=$project_dir"/data/subsidiary-files/gibbscluster-results"

sbatch_script=$project_dir"/scripts/gibbscluster-sbatch.sh"

file_list=$(find $input_dir -maxdepth 1 -type f -name "*.tsv")

for f in ${file_list[@]}; do
  f_name=$(basename $f)
  output_file=$output_dir"/output_file_"$f_name".csv.out"
  sbatch --job-name="gibbs_cl3-7-"$f_name $sbatch_script $f $output_dir $output_file
done
