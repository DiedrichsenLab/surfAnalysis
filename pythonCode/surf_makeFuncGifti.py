#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 21:32:02 2019

@author: switt
"""
import os
import numpy as np
import nibabel as nb


def makeFuncGifti(data,anatomicalStruct='CortexLeft',columnNames=[]):
    
    [N,Q] = data.shape()
    #Make columnNames if empty
    if data.size() == 0:
        for i in range(Q):
            columnNames[i] = print("col_{}".format(i))
    
    M = nb.gifti.gifti.GiftiImage()        
    
    M._meta.data[0].name = 'AnatomicalStructurePrimary'
    M._meta.data[0].value = anatomicalStruct
    M._meta.data[1].name = 'encoding'
    M._meta.data[1].value = 'XML_BASE64_GZIP'
    M._labeltable.labels[0].key = 0
    M._labeltable.labels[0].label= '???'
    M._labeltable.labels[0].red = 1.0
    M._labeltable.labels[0].green = 1.0
    M._labeltable.labels[0].blue = 1.0
    M._labeltable.labels[0].alpha = 0.0
    
    for i in range(Q):
        M.darrays[i].data = np.float32(data[:,i])
        M.darrays[i].meta.data[0].name = 'Name'
        M.darrays[i].meta.data[0].value = columnNames[i]
        M.darrays[i].dims = N
        M.darrays[i].datatype = 'NIFTI_TYPE_FLOAT32'
        M.darrays[i].Intent = 'NIFTI_INTENT_NONE'
        M.darrays[i].coordsys.dataspace = 0
        
    return(M)