#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 21:03:53 2019

@author: switt
"""

def getStruct(G):
    N = len(G._meta.data)
    anatomicalStruct = []
    for i in range(N):
        if 'AnatomicalStructurePrimary' in G._meta.data[i].name:
            anatomicalStruct = G._meta_data[i].value
            return
            