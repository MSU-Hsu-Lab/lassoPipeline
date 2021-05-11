#!/usr/bin/env bash
#SBATCH -n 1 -c 1
#SBATCH --time=143:59:00
#SBATCH --mem=512G
#SBATCH --job-name=lassoScratch
#SBATCH --output=%x-%a-%A.SLURMout
#SBATCH -a 1-5

# note hard paths have been changed.

module load GCC/8.3.0
module load Python/3.8.3
source 'PATHTOPYENV'
k=$SLURM_ARRAY_TASK_ID
echo $k

traitname=$1

OUTDIR='PARENT DIR TO OUTPUT'/$traitname/
mkdir -p $OUTDIR

genoPATH='PATH TO BEDMATRIX'

python3 lasso.pysnp.py --geno-path $genoPATH \
	--trait $traitname \
	--index-var $k \
	--output-directory $OUTDIR 
