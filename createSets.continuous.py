import numpy as np
import pandas as pd
import scipy as sp
import math as math
import time as time
import random as random
import argparse
import textwrap
import os

from pandas_plink import read_plink

def createSets(genPATH,trait,sibsPATH,whitePATH,outDIR):

    phenPATH = outDIR+trait+".pruned.txt"
    sibsPATH = sibsPATH
    whitePATH = whitePATH

    phen = pd.read_csv(phenPATH,header=0,sep=' ')
    white = pd.read_csv(whitePATH,sep=' ',header=None)
    sibs = pd.read_csv(sibsPATH,sep=' ',header=None)
    genotyped = pd.read_csv(genPATH+".fam",sep=' ',header=None)

    genLoc = phen["EID"].isin(genotyped[0])
    whiteLoc = phen["EID"].isin(white[0])
    nonWhiteLoc = ~whiteLoc & genLoc
    sibLoc = phen["EID"].isin(sibs[0]) & genLoc
    missing = phen.iloc[:,3].astype(str).eq("NA") | phen.iloc[:,3].astype(str).eq("NaN") | phen.iloc[:,3].astype(str).eq("nan")
    remove = sibLoc | nonWhiteLoc | missing | ~genLoc

    out1=outDIR+"nonWhite.txt"
    out2=outDIR+"siblings."+trait+".txt"
    phen[nonWhiteLoc].to_csv(r''+out1, sep=' ',index=False,header=False)
    phen[sibLoc].to_csv(r''+out2, sep=' ',index=False,header=False)

    subsetted = pd.DataFrame(list(range(phen.shape[0])))[-remove]
    randOrder = random.sample(range(subsetted.shape[0]),subsetted.shape[0])
    randPhen = phen.iloc[subsetted[0]].iloc[randOrder]

    N = 1000

    for i in range(1,6):
    
        print(i)
        valKeep = randPhen[(N*(i-1)):(N*i)]
        trainKeep = randPhen[~randPhen["EID"].isin(valKeep["EID"])]
    
        valSet = phen[phen["EID"].isin(valKeep["EID"])]["EID"]
        trainSet = phen[phen["EID"].isin(trainKeep["EID"])]["EID"]

        name1=outDIR+"TrainSet."+str(i)+".txt"
        name2=outDIR+"ValSet."+str(i)+".txt"
    
        trainSet.to_csv(r''+name1,sep = ' ',index=True,header=False)
        valSet.to_csv(r''+name2,sep = ' ',index=True,header=False)

    return 0
        
def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     prog='prunePheno',
                                     usage='%(prunePheno)s',
                                     description='''Prunes raw pheno to match geno.''')

    # essential arguments
    required_named = parser.add_argument_group('Required named arguments')
    required_named.add_argument('--geno-path',
                                type=str,
                                required=True,
                                help='path to genotypes')

    required_named.add_argument('--trait',
                                type=str,
                                required=True,
                                help='name of trait')

    required_named.add_argument('--sibs-path',
                                type=str,
                                required=True,
                                help='path to sibs')

    required_named.add_argument('--white-path',
                                type=str,
                                required=True,
                                help='path to whites')

    # file to
    required_named.add_argument('--output-directory',
                                type=str,
                                required=True,
                                help='Where all the output goes')

    args = parser.parse_args()
    createSets(args.geno_path,args.trait,args.sibs_path,args.white_path,args.output_directory)

exit(main())
