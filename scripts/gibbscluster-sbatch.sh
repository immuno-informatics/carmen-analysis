#!/bin/bash -l

#SBATCH --partition=partition
#SBATCH --account=account

#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --cpus-per-task=24
#SBATCH --mem=48G

#SBATCH --time=72:00:00

#SBATCH --output="/path/to/carmen-analysis/temp/slurm-logs/%x-%j.out"
#SBATCH --error="/path/to/carmen-analysis/temp/slurm-logs/%x-%j.err"

input_file=$1
output_dir=$2
output_file=$3

f_name=$(basename $input_file)

gibbscluster -f $input_file -R $output_dir -P $f_name".out-len-9" -T -j 10 -g3-7 -l9 -S1 > $output_file
