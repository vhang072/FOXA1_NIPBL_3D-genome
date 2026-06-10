# -*- coding: utf-8 -*-
"""

@author: X.Z.
"""

import numpy as np
import pandas as  pd


loops = pd.read_csv("merged.loops.txt",sep = "\t",header = 0,dtype = {'V2':np.int32,\
                                                                                                                  'V3':np.int32,\
                                                                                                                  'V5':np.int32,\
                                                                                                                  'V6':np.int32,})
bins = pd.read_csv("siCtrl.1_10000_abs.bed",sep = "\t",header = None,\
                   dtype = {1:np.int32,2:np.int32,3:np.int32})

contact = pd.read_csv("siNIPBL.1_10000.matrix",sep = "\t",\
                      header = None,dtype = {0:np.int32,1:np.int32,2:np.int32})
output = open(r"siNIPBL.1_10000.loop.counts.txt","w")    

chroms = ["chr"+str(x) for x in range(1,23)] + ["chrX","chrY"]
for i in chroms:
    print(i+"\n")
    loops_chrom = loops.loc[loops.iloc[:,0]==i.split("r")[1],]
    bins_chrom = bins.loc[bins.iloc[:,0]==i,]
    contact_chrom = contact.loc[np.array(contact.iloc[:,0]>=np.min(bins_chrom.iloc[:,3])) & np.array(contact.iloc[:,0]<=np.max(bins_chrom.iloc[:,3])) &\
                                np.array(contact.iloc[:,1]>=np.min(bins_chrom.iloc[:,3])) & np.array(contact.iloc[:,1]<=np.max(bins_chrom.iloc[:,3])),:]
    for j in range(loops_chrom.shape[0]):
        if int(np.ceil(j/100))*100 ==j:
            print("Has processed "+ str(j) + " Loops\n")
        idx1 = int(bins_chrom.loc[bins.iloc[:,1]==loops_chrom.iloc[j,1],3])
        idx2 = int(bins_chrom.loc[bins.iloc[:,1]==loops_chrom.iloc[j,4],3])
        f1 = np.array(contact_chrom.iloc[:,0]==idx1) & np.array(contact_chrom.iloc[:,1]==idx2)
        if np.sum(f1) == 1:
            contact_10k = int(contact_chrom.loc[f1,2])
        else:
            contact_10k = 0
        contact_30k = 0
        for k in range(0,3):
            for m in range(0,3):
                f1 = np.array(contact_chrom.iloc[:,0]==idx1-1+k) & np.array(contact_chrom.iloc[:,1]==idx2-1+m)
                if np.sum(f1) == 1:
                    contact_30k += int(contact_chrom.loc[f1,2])
        output.write(i + "\t" +\
                     str(loops_chrom.iloc[j,1])+ "\t" +\
                     str(loops_chrom.iloc[j,2])+ "\t" +\
                     i + "\t" +\
                     str(loops_chrom.iloc[j,4])+ "\t" +\
                     str(loops_chrom.iloc[j,5])+ "\t" +\
                     str(contact_10k)+ "\t" +\
                     str(contact_30k)+ "\n")
output.close()