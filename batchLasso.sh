#!/usr/bin/env bash
#SBATCH -n 1 -c 1
#SBATCH --time=143:59:00
#SBATCH --mem=512G
#SBATCH --job-name=lassoScratch
#SBATCH --output=%x-%a-%A.SLURMout
#SBATCH -a 1-5

module load GCC/8.3.0
module load Python/3.8.3
source /mnt/home/lellolou/myPy/python3.8/newPy/bin/activate
k=$SLURM_ARRAY_TASK_ID
echo $k

traitname=$1

OUTDIR=/mnt/home/lellolou/siblingGP/$traitname/
mkdir -p $OUTDIR

#genoPATH=/mnt/research/UKBB/hsuGroup/ukb500/genotypes/calls.merged/ukb500.calls.onlyqc
genoPATH=/mnt/home/lellolou/hsuGroup/ukb500/genotypes/calls.merged/ukb500.calls.gpsnp.biomarker3

python3 lasso.pysnp.py --geno-path $genoPATH \
	--trait $traitname \
	--index-var $k \
	--output-directory $OUTDIR 
