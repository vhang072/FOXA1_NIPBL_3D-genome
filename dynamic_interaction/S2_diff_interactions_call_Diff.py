## DESCRIPTIONS of THIS CODE (only used for intrachromosome)
# author: X.Z.
# reference Dixon et al. 2012
# (STEP 2)This code was used to perform differential interaction analysis and calculate P values, and P values of background. (See methods)
#input files:
#		1.the list of raw matrix (not be normalized) files (N*N) from chr1:chrX; two replicates of two states (4 lists)
#		2.mean_variance_number which calculated in last step. we need four files for two replicates of two states. (4 files)
#output files:
#		1.the list of P values for one foreground and two background


#note: we have mark the positions where we need to set again using "##**##".
import pandas as pd
import getopt,sys

local = 600#(6M/10k) ##**##the max distance(/bin) to calculate . This is depend on matrix resolution and sequencing depth

argumentList = sys.argv[1:]
options = "a:b:c:d:o:m:"
long_options = ["sampleA=","sampleB=","sampleC=","sampleD=","output=","maxbin="]
arguments, values = getopt.getopt(argumentList, options, long_options)
for currentArgument, currentValue in arguments:
    if currentArgument in ("-a", "--sampleA"):
        sa = currentValue
    elif currentArgument in ("-b", "--sampleB"):
        sb = currentValue
    elif currentArgument in ("-c", "--sampleC"):
        sc = currentValue
    elif currentArgument in ("-d", "--sampleD"):
        sd = currentValue
    elif currentArgument in ("-o", "--output"):
        out = currentValue
    elif currentArgument in ("-m", "--maxbin"):
        local = int(currentValue)
        
minc = 10#(mini counts) ##**## sum of count for one interaction less than this num will be ignored
import numpy as np
import math
from scipy.stats import binom
#means
file1 = np.loadtxt(sa+".mean_variance_number")##**##two replicates mean_variance_number files
lens1 = file1.shape[0]
means = np.zeros([lens1,4])
means[:,0] = file1[0:lens1,0]
file1 = np.loadtxt(sb+".mean_variance_number")
means[:,1] = file1[0:lens1,0]
file1 = np.loadtxt(sc+".mean_variance_number")
means[:,2] = file1[0:lens1,0]
file1 = np.loadtxt(sd+".mean_variance_number")
means[:,3] = file1[0:lens1,0]

#dynamic interactions
re = pd.DataFrame(index = range(int(3e8)),columns = ["chrom","bin1","bin2","sample1","sample2","sample3","sample4",\
                                                     "scaleS1","scaleS2","scaleS3","scaleS4","TrueP","BgP1","BgP2"])
ids = 0
chroms = ["chr"+str(i) for i in range(1,23)] + ["chrX"]
for i in chroms:
    data1 = np.loadtxt(sa+"_"+i+'_dense.matrix')##**## replicate1 matrix of state 1
    data2 = np.loadtxt(sb+"_"+i+'_dense.matrix')##**## replicate2 matrix of state 1
    data3 = np.loadtxt(sc+"_"+i+'_dense.matrix')##**## replicate1 matrix of state 2
    data4 = np.loadtxt(sd+"_"+i+'_dense.matrix')##**## replicate2 matrix of state 2
    lens1 = data1.shape[0]
    zrow = np.zeros([lens1,4])
    zrow[:,0] = (np.sum(data1,axis = 1)==0)
    zrow[:,1] = (np.sum(data2,axis = 1)==0)
    zrow[:,2] = (np.sum(data3,axis = 1)==0)
    zrow[:,3] = (np.sum(data4,axis = 1)==0)
    for j in range(0,lens1):
        if sum(zrow[j,0:4])==0:
            s = max(0,j-local)
            o = j
            for k in range(s,o+1):
                val = [data1[j,k],data2[j,k],data3[j,k],data4[j,k]]
                if sum(zrow[k,0:4])==0 and sum(val)>=minc:
                    p1 = sum(means[j-k,0:2])*1.0/sum(means[j-k,0:4])
                    p2 = sum(means[j-k,[0,2]])*1.0/sum(means[j-k,0:4])
                    p3 = sum(means[j-k,[0,3]])*1.0/sum(means[j-k,0:4])
                    trials = math.ceil(sum(val))
                    norm_vals = [round(x,1) for x in np.array(val,dtype = float)/means[j-k,0:4]*np.power((np.prod(means[j-k,0:4]*1.0)),(1/4))]
                    if sum(val[0:2])-1>=p1*trials:
                        pval1 = math.log10(1-binom.cdf(sum(val[0:2])-1,trials,p1)+1e-12)
                    else:
                        pval1 = -1*math.log10(binom.cdf(sum(val[0:2]),trials,p1)+1e-12)
                    if val[0]+val[2]-1>=p2*trials:
                        pval2 = math.log10(1-binom.cdf(val[0]+val[2]-1,trials,p2)+1e-12)
                    else:
                        pval2 = -1*math.log10(binom.cdf(val[0]+val[2],trials,p2)+1e-12)
                    if val[0]+val[3]-1>=p3*trials:
                        pval3 = math.log10(1-binom.cdf(val[0]+val[3]-1,trials,p3)+1e-12)
                    else:
                        pval3 = -1*math.log10(binom.cdf(val[0]+val[3],trials,p3)+1e-12)
                    re.iloc[ids,:] = [i,j,k]+val+norm_vals+[round(pval1,2),round(pval2,2),round(pval3,2)]
                    ids = ids+1

re2 = re.iloc[0:ids,:]
re2.astype({"chrom":str, "bin1":np.int32, "bin2":np.int32,"sample1":np.int32,"sample2":np.int32,\
         "sample3":np.int32,"sample4":np.int32,"scaleS1":float,"scaleS2":float,"scaleS3":float,\
         "scaleS4":float,"TrueP":float,"BgP1":float,"BgP2":float})
re2.to_csv(out+'.dynamics',sep = "\t",float_format='%.2f',index = False, header = True)

                        
                    
                    
                    
                