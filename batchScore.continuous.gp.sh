#!/usr/bin/env bash
#SBATCH -n 1 -c 10
#SBATCH --time=3:59:00
#SBATCH --mem=512G
#SBATCH --job-name=scoreCont
#SBATCH --output=%x-%a-%A.SLURMout
#SBATCH -a 1-5

module load GCC/8.3.0
module load Python/3.8.3
source /mnt/home/lellolou/myPy/python3.8/newPy/bin/activate
k=$SLURM_ARRAY_TASK_ID
echo $k

traitname=$1

OUTDIR=/mnt/home/lellolou/siblingGP/$traitname/
genoPATH=/mnt/research/UKBB/hsuGroup/ukb500/genotypes/calls.merged/ukb500.calls.gpsnp
#genoPATH=/mnt/home/lellolou/hsuGroup/ukb500/genotypes/calls.merged/ukb500.calls.gpsnp.biomarker3

#scratchPATH=/mnt/home/lellolou/scratch/ukb500.calls.gpsnp
cohortDIR=/mnt/research/UKBB/hsuGroup/ukb500/cohorts/

python3 scoreSets.continuous.py --geno-path $genoPATH \
	--trait $traitname \
	--cohort-path $cohortDIR \
	--array-id $k \
	--output-directory $OUTDIR 
