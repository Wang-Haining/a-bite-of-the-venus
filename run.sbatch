#!/bin/bash
#####  Constructed by HPC everywhere #####
#SBATCH --account=general
#SBATCH --mail-user=hw56@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --time=4-0:0:00
#SBATCH --mem=256gb
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --job-name=venus
#SBATCH --output=/N/scratch/hw56/a-bite-of-the-venus/logs/rscript-%j.out
#SBATCH --error=/N/scratch/hw56/a-bite-of-the-venus/logs/rscript-%j.err

######  Module commands #####
module load r/4.1.1

######  Job commands go below this line #####
for ((i=1; i<=456; i+=114)); do
    Rscript Execute_Terminal.R $i $((i+113))
done