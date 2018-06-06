#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Matrix of distance in phylogenetic gree for group of interested species
    Start   : 06/04/2017
    End     : 06/08/2017
'''
import os
import argparse
import shutil
from Bio import SeqIO
from ete3 import Tree
import csv
import numpy as np
import scipy.stats as st
def parser_code():

    parser = argparse.ArgumentParser(description='The purpose of this script to p')
    
    parser.add_argument("-i", "--input_tree", help="Input tree that we want to debias")
    
    parser.add_argument("-s", "--species", help="Interest species")
      
    parser.add_argument("-o", "--output", help="output csv file")                               
    return parser.parse_args()
    
def parse(handle):
    species = []
    for line in handle.readlines():
        species.append(line.strip().split(':')[0])
    return species

def get_max_mean_variance(tree):
    leaves = tree.get_leaves()
    size = len(leaves)
    current_max = 0
    current_min = 10000
    distances = []
    for i in range(size):
        for j in range(i+1,size):
            distance = leaves[i].get_distance(leaves[j])
            if distance > current_max:
                current_max = distance
            if distance < current_min:
                current_min = distance
            distances.append(distance)
    return current_min,current_max,np.mean(distances),np.var(distances)

# given tree and species, find the current min, current max, mean, and variance of all the distance in the tree
def get_distance(tree,species,current_min,current_max,mean,var):
    first_row = ['']
    first_row.extend(species)
    matrix1 =[first_row]
    matrix2 =[first_row]
    for i in range(len(species)):
        current1 = [species[i]]
        current2 = [species[i]]
        for j in range(len(species)):
            i_node = tree&species[i]
            j_node = tree&species[j]
            distance = i_node.get_distance(j_node)
            z_score = round((distance- mean)/(var**.5),4)
            normalize = round((distance-current_min)/(current_max-current_min),4)
            current1.append(str(z_score))
            current2.append(str(normalize))
        matrix1.append(current1)
        matrix2.append(current2)
    return matrix1,matrix2

          
    
if __name__ == "__main__":
    args  = parser_code()
    input_tree = args.input_tree
    output     = args.output
    species    = args.species
    species    = sorted(parse(open(species,'r')))
    tree       = Tree(input_tree,format = 1)
    current_min,current_max,mean,var = get_max_mean_variance(tree)
    matrix1,matrix2     = get_distance(tree,species,current_min,current_max,mean,var)
    
    writer=csv.writer(open(output + ".csv","w"))
    for item in matrix1:
        writer.writerow(item)
    writer.writerow([])
    writer.writerow([])
    for item in matrix2:
        writer.writerow(item)
    