# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 16:49:54 2023

@author: ABC
"""
import numpy as np
import pandas as  pd
sps = ["siCtrl.1_10000",\
       "siCtrl.2_10000",\
       "siNIPBL.1_10000",\
       "siNIPBL.2_10000"]
for sp in sps:
    bins = pd.read_csv("./"+"siCtrl.1_10000"+"_abs.bed",sep = "\t",header = None,\
                       dtype = {1:np.int32,2:np.int32,3:np.int32})
    contact = pd.read_csv("./"+sp+".matrix",sep = "\t",\
                          header = None,dtype = {0:np.int32,1:np.int32,2:np.int32})
    chroms = ["chr"+str(x) for x in range(1,23)] + ["chrX","chrY"]
    contacts_across_distances = np.zeros([2000,2],dtype=np.int32)
    contacts_across_distances[:,0] = range(0,20000000,10000)
    for i in chroms:
        print(i+"\n")
        bins_chrom = bins.loc[bins.iloc[:,0]==i,]
        contact_chrom = contact.loc[np.array(contact.iloc[:,0]>=np.min(bins_chrom.iloc[:,3])) & np.array(contact.iloc[:,0]<=np.max(bins_chrom.iloc[:,3])) &\
                                    np.array(contact.iloc[:,1]>=np.min(bins_chrom.iloc[:,3])) & np.array(contact.iloc[:,1]<=np.max(bins_chrom.iloc[:,3])),:]
        for idx in range(0,contact_chrom.shape[0]):
            if contact_chrom.iloc[idx,1] - contact_chrom.iloc[idx,0]<=1999:
                contacts_across_distances[int(abs(contact_chrom.iloc[idx,1] - contact_chrom.iloc[idx,0])),1] += contact_chrom.iloc[idx,2]
    np.savetxt("./"+sp+".totalcount.at.a.distance",\
            contacts_across_distances,fmt="%d",delimiter = "\t")   