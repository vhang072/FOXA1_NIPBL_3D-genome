# -*- coding: utf-8 -*-
"""
Created on Wed May 10 14:18:44 2023

@author: ABC
"""

import pandas as pd
import numpy as np


fdr = 0.05
minbin = 5
"""
import getopt,sys
argumentList = sys.argv[1:]
options = "i:o:m:r:"
long_options = ["input=","output=","maxbin=","res="]
arguments, values = getopt.getopt(argumentList, options, long_options)
for currentArgument, currentValue in arguments:
    if currentArgument in ("-i", "--input"):
        prefix_in = currentValue
    elif currentArgument in ("-r", "--res"):
        res = int(currentValue)
    elif currentArgument in ("-o", "--output"):
        out = currentValue
    elif currentArgument in ("-m", "--maxbin"):
        maxbin = int(currentValue)
"""
prefix_in = r"F:\Data_storage_PC\D0.loop_strength_diff\matrix_dense\20kb_ice\h0_1_20000_iced"
out = r"I:\Projects\D0.loop_strength_diff\midfiles\tad\h0_1_20000_iced"
res = 20000
maxbin = 50


chroms = ["chr" + str(x) for x in range(1,23)] + ["chrX"]

insulations = pd.DataFrame([],columns = ["chrom","start","end"]+["InsOfWin" +str(x) + "Bins" for x in range(minbin,maxbin+1) ])
boundarydelta = pd.DataFrame([],columns = ["chrom","start","end"]+["DeltaOfWin" +str(x) + "Bins" for x in range(minbin,maxbin+1) ])
boundaryout = pd.DataFrame([],columns = ["chrom","start","end","fracAcrossDiffWindows"])
tad = pd.DataFrame([],columns = ["chrom","start","end"])
for chrom in chroms:
    print(chrom + "\n")
    data = np.loadtxt(prefix_in + "_" + chrom + "_dense.matrix")
    chromlen = np.shape(data)[0]
    """
    -----------------------------------------------------------------------
    A1.calculate the true insulation scores and delta score and raw boundaries
    """
    inslen = maxbin - minbin + 1
    ins = np.zeros((chromlen, inslen))
    binsum = data.sum(axis=0)
    for i in range(0, chromlen):
        for j in range(minbin,maxbin+1):
            s1 = i-j
            o1 = i-1
            s2 = i+1
            o2 = i+j
            if s1>=0 and o2<=chromlen-1:
                f1 = binsum[s1:(o1+1)]
                f2 = binsum[s2:(o2+1)]
                if f1.sum()!=0 and f2.sum()!=0:
                    val = data[s1:(o1+1),s2:(o2+1)]
                    values2 = val[f1!=0,:]
                    values = values2[:,f2!=0]
                    ins[i,j-minbin] = values.mean()+0.1
                else:
                    ins[i,j-minbin] = np.nan
            else:
                ins[i,j-minbin] = np.nan
    for j in range(minbin,maxbin+1):
        m = ins[np.isnan(ins[:,j-minbin])==False,j-minbin].mean()
        ins[np.isnan(ins[:,j-minbin])== False,j-minbin] = np.log2(ins[np.isnan(ins[:,j-minbin])== False,j-minbin])-np.log2(m)
    #record ins,used in the following analysis
    anno = np.vstack([[np.repeat(chrom,repeats=chromlen)],[np.arange(0,chromlen*res,res)],\
              [np.arange(res,chromlen*res+res,res)]]).T
    cat_ins = pd.DataFrame(np.hstack([anno,ins]),columns = \
                           ["chrom","start","end"]+["InsOfWin" +str(x) + "Bins" for x in range(minbin,maxbin+1)])
    insulations = pd.concat([insulations,cat_ins],ignore_index = True)
    
    #delta and boundary analyses
    delta = np.zeros((chromlen,inslen))
    boundary = np.zeros((chromlen,inslen))
    boundaryM = np.zeros((chromlen,inslen))
    #calculate delta scores
    for j in range(minbin,maxbin+1):
        for i in range(0,chromlen):
            s1 = i - 3
            o1 = i + 3
            if s1>=0 and o1<= chromlen-1 and np.isnan(ins[i,j-minbin]) == False and np.isinf(ins[i,j-minbin]) == False:
                f1 = ins[s1:i,j-minbin]
                f2 = ins[(i+1):(o1+1),j-minbin]
                if np.isnan(f1).sum()+np.isinf(f1).sum()<3 and np.isnan(f2).sum()+np.isinf(f2).sum()<3:
                    val1 = f1[np.isnan(f1)==False]
                    val2 = val1[np.isinf(val1)==False]
                    val3 = f2[np.isnan(f2)==False]
                    val4 = val3[np.isinf(val3)==False]
                    delta[i,j-minbin] = val4.mean()-val2.mean()
                else:
                    delta[i,j-minbin] = np.nan
            else:
                delta[i,j-minbin] = np.nan
    #calculate boundary and boundary metrics               
    for j in range(minbin,maxbin+1):
        for i in range(0,chromlen):
            s1 = i-3
            o1 = i+3
            if np.isnan(delta[i,j-minbin]) == False and s1>=0 and o1<=chromlen-1 and \
               np.isnan(delta[i-1,j-minbin]) == False and \
               np.isnan(delta[i+1,j-minbin]) == False:
                if delta[i+1,j-minbin]>0 and delta[i-1,j-minbin]<0 and \
                   abs(delta[i,j-minbin]) <= abs(delta[i-1,j-minbin]) and \
                   abs(delta[i,j-minbin]) <= abs(delta[i+1,j-minbin]):
                    ii = i+2
                    tmp_delta = delta[i+1,j-minbin]
                    if np.isnan(delta[ii,j-minbin]) == False:
                        while delta[ii,j-minbin] > tmp_delta:
                            tmp_delta = delta[ii,j-minbin]
                            ii +=1
                            if ii > chromlen-1:
                                break
                            if ii <= chromlen -1:
                                if np.isnan(delta[ii,j-minbin]) == True:
                                    break
                    ii = i-2
                    tmp_delta2 = delta[i-1,j-minbin]
                    if np.isnan(delta[ii,j-minbin]) == False:
                        while delta[ii,j-minbin] < tmp_delta2:
                            tmp_delta2 = delta[ii,j-minbin]
                            ii -=1
                            if ii < 0:
                                break
                            if ii >=0:
                                if np.isnan(delta[ii,j-minbin]) == True:
                                    break
                    boundary[i,j-minbin] = 1
                    boundaryM[i,j-minbin] = tmp_delta -tmp_delta2
    """
    ------------------------------------------------------------------
    A2.build a shuffled matrix and calculate the background insulations
    
    """
    sign = np.zeros((chromlen,chromlen))
    sign[binsum==0,:] = 1
    sign[:,binsum ==0] = 1
    bg = np.zeros((chromlen,chromlen))
    for i in range(0,maxbin*2+1):
        signd= np.diag(sign,k=i)
        orignd = np.diag(data,k=i)
        nowd = orignd[signd == 0]
        np.random.shuffle(nowd)
        tmp = np.zeros((chromlen-i))
        tmp[signd==0] = nowd
        if i==0:
            bg += np.diag(tmp,k=i)
        else:
            bg += np.diag(tmp,k=i)
            bg += np.diag(tmp,k=-1*i)
    ins_bg = np.zeros((chromlen,inslen))
    for i in range(0,chromlen):
        for j in range(minbin,maxbin+1):
            s1 = i-j
            o1 = i-1
            s2 = i+1
            o2 = i+j
            if s1>=0 and o2<=chromlen-1:
                f1 = binsum[s1:(o1+1)]
                f2 = binsum[s2:(o2+1)]
                if f1.sum()!=0 and f2.sum()!=0:
                    val = bg[s1:(o1+1),s2:(o2+1)]
                    values2 = val[f1!=0,:]
                    values = values2[:,f2!=0]
                    ins_bg[i,j-minbin] = values.mean()+0.1
                else:
                    ins_bg[i,j-minbin] = np.nan
            else:
                ins_bg[i,j-minbin] = np.nan
    for j in range(minbin,maxbin+1):
        m = ins_bg[np.isnan(ins_bg[:,j-minbin])==False,j-minbin].mean()
        ins_bg[np.isnan(ins_bg[:,j-minbin])== False,j-minbin] = np.log2(ins_bg[np.isnan(ins_bg[:,j-minbin])== False,j-minbin])-np.log2(m)            
        #record ins_bg,used in the following analysis

    #delta and boundary analyses
    delta_bg = np.zeros((chromlen,inslen))
    boundary_bg = np.zeros((chromlen,inslen))
    boundaryM_bg = np.zeros((chromlen,inslen))
    #calculate delta scores
    for j in range(minbin,maxbin+1):
        for i in range(0,chromlen):
            s1 = i - 3
            o1 = i + 3
            if s1>=0 and o1<= chromlen-1 and np.isnan(ins_bg[i,j-minbin]) == False and np.isinf(ins_bg[i,j-minbin]) == False:
                f1 = ins_bg[s1:i,j-minbin]
                f2 = ins_bg[(i+1):(o1+1),j-minbin]
                if np.isnan(f1).sum()+np.isinf(f1).sum()<3 and np.isnan(f2).sum()+np.isinf(f2).sum()<3:
                    val1 = f1[np.isnan(f1)==False]
                    val2 = val1[np.isinf(val1)==False]
                    val3 = f2[np.isnan(f2)==False]
                    val4 = val3[np.isinf(val3)==False]
                    delta_bg[i,j-minbin] = val4.mean()-val2.mean()
                else:
                    delta_bg[i,j-minbin] = np.nan
            else:
                delta_bg[i,j-minbin] = np.nan
    #calculate boundary and boundary metrics               
    for j in range(minbin,maxbin+1):
        for i in range(0,chromlen):
            s1 = i-3
            o1 = i+3
            if np.isnan(delta_bg[i,j-minbin]) == False and s1>=0 and o1<=chromlen-1 and \
               np.isnan(delta_bg[i-1,j-minbin]) == False and \
               np.isnan(delta_bg[i+1,j-minbin]) == False:
                if delta_bg[i+1,j-minbin]>0 and delta_bg[i-1,j-minbin]<0 and \
                   abs(delta_bg[i,j-minbin]) <= abs(delta_bg[i-1,j-minbin]) and \
                   abs(delta_bg[i,j-minbin]) <= abs(delta_bg[i+1,j-minbin]):
                    ii = i+2
                    tmp_delta = delta_bg[i+1,j-minbin]
                    if np.isnan(delta_bg[ii,j-minbin]) == False:
                        while delta_bg[ii,j-minbin] > tmp_delta:
                            tmp_delta = delta_bg[ii,j-minbin]
                            ii +=1
                            if ii > chromlen-1:
                                break
                            if ii <= chromlen -1:
                                if np.isnan(delta_bg[ii,j-minbin]) == True:
                                    break
                    ii = i-2
                    tmp_delta2 = delta_bg[i-1,j-minbin]
                    if np.isnan(delta_bg[ii,j-minbin]) == False:
                        while delta_bg[ii,j-minbin] < tmp_delta2:
                            tmp_delta2 = delta_bg[ii,j-minbin]
                            ii -=1
                            if ii < 0:
                                break
                            if ii >=0:
                                if np.isnan(delta_bg[ii,j-minbin]) == True:
                                    break
                    boundary_bg[i,j-minbin] = 1
                    boundaryM_bg[i,j-minbin] = tmp_delta -tmp_delta2
    """
    ------------------------------------------------------------------
    A3.control FDR and identify ture boundaries
    
    """
    boundary_FDR = np.zeros((chromlen,inslen))
    boundaryM_FDR = np.zeros((chromlen,inslen))
    for j in range (minbin,maxbin+1):
        p1 = boundaryM[boundary[:,j-minbin]==1,j-minbin]
        p2 = boundaryM_bg[boundary_bg[:,j-minbin]==1,j-minbin]
        p = np.append(p1,p2)
        lens = np.shape(p1)[0]+np.shape(p2)[0]
        index = np.zeros((lens))
        index[0:(np.shape(p1)[0])] = 1
        p_index = np.argsort(p)
        p_sort = p[p_index]
        index_sort = index[p_index]
        fdr_real = 1
        i = 0
        while fdr_real>fdr:
            i +=1
            fdr_real = 1.0*((index_sort[i:]==0).sum())/(lens-i)
            if i >len(p2)+len(p1)/3*2:
                print("FDR could not meet the requirment; " + "FDR_real =" +str(fdr_real)+"\n")
                break
        threhold = p_sort[i]
        boundary_FDR[:,j-minbin] = boundaryM[:,j-minbin]>=threhold
        boundaryM_FDR[boundaryM[:,j-minbin]>=threhold,j-minbin] = boundaryM[boundaryM[:,j-minbin]>=threhold,j-minbin]
        
    cat_delta = pd.DataFrame(np.hstack([anno,boundaryM_FDR]),columns = \
                           ["chrom","start","end"]+["DeltaOfWin" +str(x) + "Bins" for x in range(minbin,maxbin+1)])
    boundarydelta = pd.concat([boundarydelta,cat_delta],ignore_index = True)
    """
    ------------------------------------------------------------------
    A4.merge possible boundaries across different parameters
     
    """       
    boundary_num = boundary_FDR.sum(axis = 1) + np.random.random(chromlen)
    boundary_merged = np.zeros((chromlen,1))
    for i in range(3,chromlen-3):
        if boundary_num[i] > boundary_num[i-3] and \
           boundary_num[i] > boundary_num[i-2] and \
           boundary_num[i] > boundary_num[i-1] and \
           boundary_num[i] > boundary_num[i+3] and \
           boundary_num[i] > boundary_num[i+2] and \
           boundary_num[i] > boundary_num[i+1] and \
           sum(np.floor(boundary_num[(i-3):(i+4)])) > inslen/5:
            boundary_merged[i,0] = sum(np.floor(boundary_num[(i-3):(i+4)]))/inslen
    final_boundaries = np.zeros((sum(boundary_merged[:,0]!=0),4),dtype = "<U11")
    print("final boundaries size:" + str(final_boundaries.shape))
    idx = 0
    for i in range(0,chromlen):
        if boundary_merged[i,0] !=0:
            final_boundaries[idx,:] = np.array([chrom,i*res,i*res+res,boundary_merged[i,0]])
            idx += 1
    boundary_chrom = pd.DataFrame(final_boundaries,columns = ["chrom","start","end","fracAcrossDiffWindows"])
    boundaryout = pd.concat([boundaryout,boundary_chrom],ignore_index = True)
    """
    ------------------------------------------------------------------
    A5.derive TADs according to boundary
     
    """      
    pos1 = (np.where(boundary_merged[:,0]!=0))[0]
    flag = 0
    for i in range(0,pos1.shape[0]-1):
        s = pos1[i]
        o = pos1[i+1]
        if (binsum[s+1:o]!=0).sum()>=(o-s-1)*0.8 and flag==0:
            tad_chrom = np.array([chrom,(s+1)*res,o*res])
            flag = 1
        elif (binsum[s+1:o]!=0).sum()>=(o-s-1)*0.8 and flag==1:
            tad_chrom = np.vstack((tad_chrom,np.array([chrom,(s+1)*res,o*res])))
    tad_chrom_df = pd.DataFrame(tad_chrom,columns = ["chrom","start","end"])
    tad = pd.concat([tad,tad_chrom_df],ignore_index = True)
tad.to_csv(out+".tad.txt",sep = "\t",header = True, index = False)
insulations.to_csv(out+".insulation.score.txt",sep = "\t",header = True, index = False)
boundarydelta.to_csv(out+".insulation.delta.txt",sep = "\t",header = True, index = False)
boundaryout.to_csv(out+".boundary.txt",sep = "\t",header = True, index = False)
