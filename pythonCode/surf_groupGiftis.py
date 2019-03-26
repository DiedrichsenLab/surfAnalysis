#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 08:36:34 2019

@author: switt
"""
import os
import numpy as np
import nibabel as nb
import surf_makeFuncGifti
import surf_getGiftiColumnNames
import surf_getGiftiAnatomicalStruct


def surf_groupGiftis(filelist,inputcol=np.array(),outcolnames=[],outfilenames=[],outfilenamePattern='%s.func.gii', groupsummary=[],groupstats=???):
    
    replaceNans = 0
    
    [a,b,ext] = os.path.split(P[1,:])
    
    with open(filelist, 'r') as f:
        P = f.readlines()
    
    for i in range(len(P,1)):
        filename = (P[i,:]).strip()
        if os.path.isfile(filename)=False:
            print("File {} does not exist; replacing with NaNs".format(filename))
            numCols[i]=np.nan
            numRows[i]=np.nan
        else:
            M[i]=nb.load(filename)
            da = M[i].darrays
            [numRows[i],numCols[i]]=[int(da[0].data.shape[0]),len(da))]
            
    #Check if metric files are same format
    if np.var(numCols[~np.isnan(numCols)])>0:
        print("Number of columns is not the same")
    if np.var(numRows[~np.isnan(numRows)])>0:
        print("Number of vertices is not the same")
        
    #Determine the number of columns
    if inputcol.size == 0:
        inputcol = np.arange(np.amin(numCols))
        
    numrows = np.nanmean(numRows)
    
    #Get the input column names
    inputColumnNames = surf_getGiftiColumnNames.getColumnNames(M[0])
    
    #Get main anatomical structure
    anatStruct = surf_getGiftiAnatomicalStruct.getStruct(M[0])
    
    #Determine output column names
    if outcolnames.size == 0:
        for i in range(inputcol):
            [a,b,c] = os.path.split(P[col,:])
            outcolnames[i] = b
            
    #Reorder into new metric files
    for file in range(inputcol):
        OutData = np.empty((numrows,len(M))*np.nan
        for col in range(len(M)):
            if len(M[col]) ~= 0:
                Outdata[:,col] = M[col].darrays(:,inputcol[file])
        if (len(outfilenames) == 0 or len(outfilenames)<file):
            outfilenames[file] = print("{}.func.gii".format(inputColumnNames[inputcol[file]]))
        O = surf_makeFuncGifti.makeFuncGifti(OutData, columnNames=outcolnames, anattomicalStruct=anatStruct)
        nb.save(O,outfilenames[file])
        SumData[:,file] = groupstats(OutData)
        
    #If Summary file name is given
    if len(groupsummary) ~= 0:
        O = surf_makeFuncGifti(SumData,columNames=outcolnames,anatomicalStruct=anatStruct)
        nb.save(O,groupsummary)
        
    
