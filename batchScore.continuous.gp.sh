#!/usr/bin/env bash
#SBATCH -n 1 -c 10
#SBATCH --time=3:59:00
#SBATCH --mem=512G
#SBATCH --job-name=scoreCont
#SBATCH --output=%x-%a-%A.SLURMout
#SBATCH -a 1-5

module load GCC/8.3.0
module load Python/3.8.3
source 'PATH TO PY ENV'
k=$SLURM_ARRAY_TASK_ID
echo $k

traitname=$1

OUTDIR='PARENT OUT DIR'/$traitname/
genoPATH='PATH TO BED MATRIX'

cohortDIR='PATH WHERE COHORTS ARE DEFINED'

python3 scoreSets.continuous.py --geno-path $genoPATH \
	--trait $traitname \
	--cohort-path $cohortDIR \
	--array-id $k \
	--output-directory $OUTDIR 
