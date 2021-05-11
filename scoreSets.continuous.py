import numpy as np
import pandas as pd
import scipy as sp
import math as math
import time as time
import argparse
import textwrap
import os
from pandas_plink import read_plink
from sklearn import linear_model
import sklearn as skl
from pysnptools.snpreader import Bed
            
def scoreSets(genopath,trait,cohortpath,arrayid,outputdirectory):
    
    index = arrayid
    genPATH = genopath
    phenPATH = outputdirectory+trait+".pruned.txt"
    gwasPATH = outputdirectory+"gwas."+str(trait)+"."+str(index)+".csv"
    trainPATH = outputdirectory+"TrainSet."+str(index)+".txt"
    lamPATH = outputdirectory+"lasso.lambdas."+str(trait)+"."+str(index)+".txt"
    betaPATH = outputdirectory+"lasso.betas."+str(trait)+"."+str(index)+".txt"
    corrPATH = outputdirectory+"cor."+str(trait)+"."+str(index)+".txt"
    print(outputdirectory)
    print(phenPATH)
    
    valPATH = outputdirectory+"ValSet."+str(index)+".txt"
    sibPATH = outputdirectory+"siblings."+str(trait)+".txt"
    testPATH = outputdirectory+"nonWhite.txt"
    whitePATH = cohortpath+"white.report.txt"
    asianPATH = cohortpath+"asian.report.txt"
    chinesePATH = cohortpath+"chinese.report.txt"
    blackPATH = cohortpath+"black.report.txt"

    scoreValOUT = outputdirectory+"val.score."+str(trait)+"."+str(index)+".txt"
    scoreSibsOUT = outputdirectory+"sib.score."+str(trait)+"."+str(index)+".txt"
    scoreAsianOUT = outputdirectory+"asian.score."+str(trait)+"."+str(index)+".txt"
    scoreChineseOUT = outputdirectory+"chinese.score."+str(trait)+"."+str(index)+".txt"
    scoreBlackOUT = outputdirectory+"black.score."+str(trait)+"."+str(index)+".txt"
    pruneBetaOUT = outputdirectory+"beta."+str(trait)+"."+str(index)+".txt"
    
#    (bim,fam,G) = read_plink(genPATH)
#    fam = pd.DataFrame(fam.values)
#    bim = pd.DataFrame(bim.values)

    G = Bed(genPATH,count_A1=False)
    fam = pd.read_csv(genPATH+".fam",header=None,sep=' ')
    bim = pd.read_csv(genPATH+".bim",header=None,sep='\t')

    phen = pd.read_csv(phenPATH,header=0,sep=' ')
    beta = pd.read_csv(betaPATH,header=None,sep=' ')

    # ID missing phenotypes
    keep = ~phen.iloc[:,3].isnull()
    
    names = beta[1]
    subsetP = bim[1].isin(names)
    subsetP = np.stack(pd.DataFrame(list(range(bim.shape[0])))[subsetP].values,axis=1)[0]

    print("read sets...")
    
    valSet = pd.read_csv(valPATH,header=None,sep=' ')
    sibs = pd.read_csv(sibPATH,header=None,sep=' ')
    nonwhite = pd.read_csv(testPATH,header=None,sep=' ')
    white = pd.read_csv(whitePATH,header=None,sep=' ')
    asian = pd.read_csv(asianPATH,header=None,sep=' ')
    chinese = pd.read_csv(chinesePATH,header=None,sep=' ')
    black = pd.read_csv(blackPATH,header=None,sep=' ')

    valN = valSet[0].values

    nonwhiteKeep = fam[0].astype(int).isin(nonwhite[0]) #& keep
    whiteKeep = fam[0].astype(int).isin(white[0]) & keep
    asianKeep = fam[0].astype(int).isin(asian[0]) #& keep
    chineseKeep = fam[0].astype(int).isin(chinese[0]) #& keep
    blackKeep = fam[0].astype(int).isin(black[0]) #& keep
    sibsKeep = fam[0].astype(int).isin(sibs[0]) & keep

    print("Id indices")
    
    sibsN = np.stack(pd.DataFrame(list(range(fam.shape[0])))[sibsKeep & whiteKeep].values,axis=1)[0]
    asianN = np.stack(pd.DataFrame(list(range(fam.shape[0])))[nonwhiteKeep & asianKeep].values,axis=1)[0]
    chineseN = np.stack(pd.DataFrame(list(range(fam.shape[0])))[nonwhiteKeep & chineseKeep].values,axis=1)[0]
    blackN = np.stack(pd.DataFrame(list(range(fam.shape[0])))[nonwhiteKeep & blackKeep].values,axis=1)[0]

    print("Calling into memory...",flush=True)
    # this will subset the bed matrix
    # and then actually load it into memory (.compute())
    # find missing values in the first column and count how many
    t = time.time()
    subG = G[:,subsetP].read().val
    valG = subG[valN,:]
    sibG = subG[sibsN,:]
    asianG = subG[asianN,:]
    chineseG = subG[chineseN,:]
    blackG = subG[blackN,:]
    subG = 0.0
    elapsed = time.time() - t
    print(elapsed)
    # note, 1000 columns loads basically as quick as 1
    print("Final shapes:")
    print(valG.shape)
    print(sibG.shape)
    print(asianG.shape)
    print(chineseG.shape)
    print(blackG.shape)

    print("Calc means")
    # calculate column means with no missing values
    # nanmean calculates mean skipping nan
    center = beta[6]
    
    print("NA repl")
    # na replacement
    missing = np.argwhere(np.isnan(valG))
    for row in range(0,missing.shape[0]):
        ind1 = missing[row,0]
        ind2 = missing[row,1]
        valG[ind1,ind2] = center[ind2]
        
    missing = np.argwhere(np.isnan(sibG))
    for row in range(0,missing.shape[0]):
        ind1 = missing[row,0]
        ind2 = missing[row,1]
        sibG[ind1,ind2] = center[ind2]

    missing = np.argwhere(np.isnan(asianG))
    for row in range(0,missing.shape[0]):
        ind1 = missing[row,0]
        ind2 = missing[row,1]
        asianG[ind1,ind2] = center[ind2]

    missing = np.argwhere(np.isnan(chineseG))
    for row in range(0,missing.shape[0]):
        ind1 = missing[row,0]
        ind2 = missing[row,1]
        chineseG[ind1,ind2] = center[ind2]

    missing = np.argwhere(np.isnan(blackG))
    for row in range(0,missing.shape[0]):
        ind1 = missing[row,0]
        ind2 = missing[row,1]
        blackG[ind1,ind2] = center[ind2]

        
    print("Center")
    for col in range(0,valG.shape[1]):
        valG[:,col] = valG[:,col] - center[col]
        sibG[:,col] = sibG[:,col] - center[col]
        asianG[:,col] = asianG[:,col] - center[col]
        chineseG[:,col] = chineseG[:,col] - center[col]
        blackG[:,col] = blackG[:,col] - center[col]

    valScore = np.dot(valG,(beta.iloc[:,7:beta.shape[1]]))
    valPhen = phen.iloc[valN,3]
    sibScore = np.dot(sibG,(beta.iloc[:,7:beta.shape[1]]))
    sibPhen = phen.iloc[sibsN,3]
    asianScore = np.dot(asianG,(beta.iloc[:,7:beta.shape[1]]))
    asianPhen = phen.iloc[asianN,3]
    chineseScore = np.dot(chineseG,(beta.iloc[:,7:beta.shape[1]]))
    chinesePhen = phen.iloc[chineseN,3]
    blackScore = np.dot(blackG,(beta.iloc[:,7:beta.shape[1]]))
    blackPhen = phen.iloc[blackN,3]

    aucVal = np.zeros(valScore.shape[1])

    print("SIZE VALSCORE:")
    print(valScore.shape)
    
    for k in range(1,valScore.shape[1]):
#        aucVal[k] = skl.metrics.roc_auc_score(valPhen, valScore[:,k])
        aucVal[k] = np.corrcoef(valPhen,valScore[:,k])[0,1]

    pd.DataFrame(aucVal).to_csv(corrPATH,header=False,index=False)
        
    best = np.argmax(aucVal)
    print("Val R ("+str(trait)+"): "+str(aucVal[best]))
    print("Sibling (white) R: "+str(np.corrcoef(sibPhen,sibScore[:,best])[0,1]))
#    print("Asian R: "+str(np.corrcoef(asianPhen,asianScore[:,best])[0,1]))
#    print("Black R: "+str(np.corrcoef(blackPhen,blackScore[:,best])[0,1]))
#    print("Chinese R: "+str(np.corrcoef(chinesePhen,chineseScore[:,best])[0,1]))

    keepcols = [0,1,2,3,4,5,6,7+best]
    outbeta = beta.iloc[:,keepcols]
    outbeta.to_csv(r''+pruneBetaOUT,sep=' ',index=False,header=False)
    
    metadat = phen.iloc[sibsN,0:4]
    metadat = metadat.reset_index(drop=True)
    finalScore = pd.Series(sibScore[:,best])
    out = pd.concat([metadat,finalScore],ignore_index=True,axis=1)
    out.to_csv(r''+scoreSibsOUT,sep=' ',index=False,header=False)

    metadat = phen.iloc[valN,0:4]
    metadat = metadat.reset_index(drop=True)
    finalScore = pd.Series(valScore[:,best])
    out = pd.concat([metadat,finalScore],ignore_index=True,axis=1)
    out.to_csv(r''+scoreValOUT,sep=' ',index=False,header=False)

    metadat = phen.iloc[asianN,0:4]
    metadat = metadat.reset_index(drop=True)
    finalScore = pd.Series(asianScore[:,best])
    out = pd.concat([metadat,finalScore],ignore_index=True,axis=1)
    out.to_csv(r''+scoreAsianOUT,sep=' ',index=False,header=False)

    metadat = phen.iloc[chineseN,0:4]
    metadat = metadat.reset_index(drop=True)
    finalScore = pd.Series(chineseScore[:,best])
    out = pd.concat([metadat,finalScore],ignore_index=True,axis=1)
    out.to_csv(r''+scoreChineseOUT,sep=' ',index=False,header=False)

    metadat = phen.iloc[blackN,0:4]
    metadat = metadat.reset_index(drop=True)
    finalScore = pd.Series(blackScore[:,best])
    out = pd.concat([metadat,finalScore],ignore_index=True,axis=1)
    out.to_csv(r''+scoreBlackOUT,sep=' ',index=False,header=False)
    
    return 0

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     prog='scoreSets',
                                     usage='%(scoreSets)s',
                                     description='''Choosing lambda and scores.''')

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

    required_named.add_argument('--cohort-path',
                                type=str,
                                required=True,
                                help='path to whites')

    required_named.add_argument('--array-id',
                                type=str,
                                required=True,
                                help='path to whites')

    # file to
    required_named.add_argument('--output-directory',
                                type=str,
                                required=True,
                                help='Where all the output goes')

    args = parser.parse_args()
    scoreSets(args.geno_path,args.trait,args.cohort_path,args.array_id,args.output_directory)

exit(main())
