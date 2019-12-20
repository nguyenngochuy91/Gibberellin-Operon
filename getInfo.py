#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 12:55:20 2019

@author: huyn
"""

def retrieveInfo(textFile):
    handle = open(textFile,"r")
    lines = handle.readline()
    d = {}
    for line in lines:
        line = line.strip().split(":")
        name = line[0]
        geneBlock = line[1]
        d[name]   = []
        if "a" in geneBlock:
            d[name].append("a")
        if "k" in geneBlock:
            d[name].append("k")
        d[name] = ",".join(d[name])

def parseFile(infile):
    handle = open(infile,"r")
    