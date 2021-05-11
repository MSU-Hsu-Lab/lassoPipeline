import numpy as np
import pandas as pd
import scipy as sp
import math as math
import time as time
import argparse
import textwrap
import os

from pandas_plink import read_plink

def prunePheno(genPATH,phenDIR,trait,outDIR):

    outPATH = outDIR+trait+".pruned.txt"
    phenPATH = phenDIR+trait+".txt"

    (bim,fam,G) = read_plink(genPATH)
    fam = pd.DataFrame(fam.values)
    bim = pd.DataFrame(bim.values)

    phen = pd.read_csv(phenPATH,header=0,sep=' ')
    phen = phen[phen["EID"].isin(fam[0])]
    phen = phen.reset_index(drop=True)
    
    phen.to_csv(r''+outPATH,sep=' ',index=False,header=True,na_rep='NA')
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

    required_named.add_argument('--pheno-path',
                                type=str,
                                required=True,
                                help='path to phenotypes')

    required_named.add_argument('--trait',
                                type=str,
                                required=True,
                                help='name of trait')


    # file to
    required_named.add_argument('--output-directory',
                                type=str,
                                required=True,
                                help='Where all the output goes')

    args = parser.parse_args()

    prunePheno(args.geno_path,args.pheno_path,args.trait,args.output_directory)

exit(main())
