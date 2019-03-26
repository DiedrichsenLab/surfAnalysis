#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 21:11:31 2019

@author: switt
"""

def getColumnNames(G):
    N = len(G.darrays)
    for n in range(N):
        for i in range(len(G.darrays[n].meta.data)):
            if 'Name' in G.darrays[n].meta.data[i].name:
                names[n] = G.darrays[n].meta.data[i].value
                break