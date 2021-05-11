# run example
# bash runPipeline.continous.gp.sh bmi.genWhite

source /mnt/home/lellolou/myPy/bin/activate

traitname=$1
echo $traitname

OUTDIR=/mnt/home/lellolou/siblingGP/$traitname/
PLINK=/mnt/home/lellolou/programs/plink1.9/plink
mkdir -p $OUTDIR

# adjust the phenoDIR / genoDIR
phenoDIR=/mnt/research/UKBB/hsuGroup/ukb500/phenotypes/bioAssays/phenFiles/
#genoPATH=/mnt/research/UKBB/hsuGroup/ukb500/genotypes/calls.merged/ukb500.calls.gpsnp
genoPATH=/mnt/research/UKBB/hsuGroup/ukb500/genotypes/calls.merged/ukb500.calls.gpsnp.biomarker3

cohortDIR=/mnt/research/UKBB/hsuGroup/ukb500/cohorts/
whitePATH=$cohortDIR"white.report.txt"
sibsPATH=$cohortDIR"eid.withinSibPairs.txt"

echo Pruning..
python3 prunePheno.py --geno-path $genoPATH \
       --pheno-path $phenoDIR \
       --trait $traitname \
       --output-directory $OUTDIR 

echo Subsetting...
python3 createSets.continuous.py --geno-path $genoPATH \
       --trait $traitname \
       --sibs-path $sibsPATH \
       --white-path $whitePATH \
       --output-directory $OUTDIR 

PHENBASE=$OUTDIR/$traitname.pruned.txt
cd $OUTDIR
awk '{print $1, $1, $4}' $PHENBASE | tail -n +2 > $traitname.phen.txt

echo GWAS...
for k in {1..5}
do
    TRAINSET=$OUTDIR/TrainSet.$k.txt
    awk '{print $2, $2}' $TRAINSET > train.$k.txt
    $PLINK --bfile $genoPATH --keep train.$k.txt --pheno $traitname.phen.txt --assoc --out gwas.$traitname.$k
    cat gwas.$traitname.$k.qassoc | tr -s '[:blank:]' ',' > gwas.$traitname.$k.csv
    rm train.$k.txt
done

#sbatch batchLasso.sh $traitname

#for k in 1
#do
#    python3 lasso.py --geno-path $genoPATH \
#       --trait $traitname \
#       --index-var $k \
#       --output-directory $OUTDIR 
#done
