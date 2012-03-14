#!/usr/bin/env python
# encoding: utf-8
"""
prefutil.py

Created by Eric on 2009-12-28.
Copyright (c) 2009 Eric Brochu. All rights reserved.

Helper functions for the preferences.
"""
from __future__ import division

import pdb

from numpy import argmax


def query2prefs(query, f, bestDegree=0):
    """
    Given a sequence of query points and a function to evaluate, return a set 
    of preferences.
    """
    try:
        # for now, find the biggest gap and use that as the preference point
        D = zip([f(q) for q in query], query)
        D.sort(key=lambda x:x[0])
    
        bpoint = argmax([D[i+1][0]-D[i][0] for i in xrange(len(query)-1)]) + 1
    
        # everything on the right of the break point is preferred to 
        # everything on the left
        prefs = []
        for i, (yu, xu) in enumerate(D[:bpoint]):
            for j, (yv, xv) in enumerate(D[bpoint:]):
                if i==0 and j==len(D)-bpoint-1:
                    prefs.append((xv, xu, bestDegree))
                else:
                    prefs.append((xv, xu, 0))
    except e:
        print e
        pdb.set_trace()
    return prefs
    


    