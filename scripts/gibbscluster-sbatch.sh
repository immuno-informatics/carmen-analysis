#!/bin/bash -l

#SBATCH --partition=partition
#SBATCH --account=account

#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH --cpus-per-task=24
#SBATCH --mem=48G

#SBATCH --time=72:00:00

#SBATCH --output="/path/to/the/main/project/dir/temp/slurm-logs/%x-%j.out"
#SBATCH --error="/path/to/the/main/project/dir/temp/slurm-logs/%x-%j.err"

input_file=$1

f_name=$(basename $input_file)

gibbscluster -f $input_file -P $f_name".out-len-9" -T -j 10 -g3-7 -l9 -S1 > "output_file_"$f_name".csv.out"
