#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Given the distance of nodes of rpoB and operon tree, provide acomparison
    Start   : 06/04/2017
    End     : 06/08/2017
'''

import os
import argparse
import csv
import numpy as np
import scipy.stats as st
import time
def parser_code():

    parser = argparse.ArgumentParser(description='Given the distance of nodes between 2 target tree, provide a comparison')
    
    parser.add_argument("-t1", "--tree1", help="Species tree matrix distance(rpoB.csv)")
    
    parser.add_argument("-t2", "--tree2", help="Operon or gene tree matrix distance(operon_distance.csv)")
    
    parser.add_argument("-t","--threshold",help="threshold to be treated as HGT (default as .3)",default = .3)
      
    parser.add_argument("-o", "--output", help="output csv file")                               
    return parser.parse_args()
# given the matrix of distance of a given tree, parse into 2 list of zscore and norm score
def parse(handle):
    z_score = []
    norm_score =[]
    lines = handle.readlines()
    column_size = len(lines[0].split(','))

    for i in range(len(lines)):
        if len(lines[i])<column_size:
            i+=2
            break
        line = lines[i].strip("\r\n")
        z_score.append(line.split(','))
    for j in range(i,len(lines)):
        line = lines[j].strip("\r\n")
        norm_score.append(line.split(','))        
    return z_score,norm_score

'''
    function:given the z score and norm score for rpoB tree and operon tree, provide a comparison 
    For z score, provide bit score and percent difference, bit score is 1 if rpob distance is higher
    For norm score, provide bit score and percent difference, bit score is 1 if rpob distance is higher
    input:4 lists
    output: 4 lists
''' 
def compare(rpoB_z,rpoB_norm,operon_norm,operon_z,threshold):
    compare_z_score_bit = [[]]
    compare_norm_score_bit = [[]]
    compare_z_score = [[]]
    compare_norm_score = [[]]
    HGT = set()
    names = rpoB_z[0]
    for item in rpoB_z[0]:
        compare_z_score_bit[0].append(item)
        compare_norm_score_bit[0].append(item)
        compare_z_score[0].append(item)
        compare_norm_score[0].append(item)
    for i in range(1,15):
        rpoB_z_scores = rpoB_z[i]
        rpoB_norm_scores = rpoB_norm[i]
        operon_norm_scores = operon_norm[i]
        operon_z_scores = operon_z[i]
        current_name = rpoB_z_scores[0]
        current_z =[current_name]
        current_norm = [current_name]
        current_z_bit = [current_name]
        current_norm_bit = [current_name]
        for j in range(1,15):
            # calculate area difference z score
            rpoB_z_score =float(rpoB_z_scores[j])
            operon_z_score = float(operon_z_scores[j])
            difference_z = round(st.norm.cdf(rpoB_z_score)-st.norm.cdf(operon_z_score),4)
            if difference_z >threshold:
                name1 = names[i]
                name2 = names[j]
                if (name1,name2,str(difference_z)) not in HGT and (name2,name1,str(difference_z)) not in HGT:
                    HGT.add((name1,name2,str(difference_z)))
            if difference_z > 0: 
                difference_z_bit = 1  
            else:
                difference_z_bit = 0
            current_z_bit.append(str(difference_z_bit))
            current_z.append(str(difference_z))
            # calculate difference in norm score
            rpoB_norm_score = float(rpoB_norm_scores[j])
            operon_norm_score = float(operon_norm_scores[j])
            difference_norm = rpoB_norm_score- operon_norm_score
            if difference_norm >0:
                difference_norm_bit= 1
            else:
                difference_norm_bit = 0
            current_norm.append(str(difference_norm))
            current_norm_bit.append(str(difference_norm_bit))

        compare_z_score.append(current_z)
        compare_z_score_bit.append(current_z_bit)
        compare_norm_score.append(current_norm)
        compare_norm_score_bit.append(current_norm_bit)  

    return compare_z_score,compare_norm_score,compare_z_score_bit,compare_norm_score_bit,HGT

# write protential HGT
def write_HGT(threshold,HGT):
    writer=csv.writer(open("HGT_{}.csv".format(threshold),"w"))
    dic = {}
    # parse the accesion file
    infile = open("accession_to_common.txt","r")
    for line in infile.readlines():
        line = line.strip().split(',')
        dic[line[0]]= line[0]+','+line[1]
    for trios in sorted(HGT, key = lambda x : x[2]):
        writer.writerow([dic[trios[0]],dic[trios[1]],trios[2]])
# write from list into csv
def writing(filename,matrix1,matrix2):
    writer=csv.writer(open(filename + ".csv","w"))
    for item in matrix1:
        writer.writerow(item)
    writer.writerow([])
    writer.writerow([])
    for item in matrix2:
       writer.writerow(item)
       
if __name__ == "__main__":
    args = parser_code()
    rpoB_tree = args.tree1
    operon_tree = args.tree2
    output  = args.output
    threshold       = float(args.threshold)
    rpoB_z,rpoB_norm = parse(open(rpoB_tree,'r'))
    operon_z,operon_norm = parse(open(operon_tree,'r'))
    compare_z_score,compare_norm_score,compare_z_score_bit,compare_norm_score_bit,HGT = compare(rpoB_z,rpoB_norm,operon_norm,operon_z,threshold)
    write_HGT(threshold,HGT)    
    writing(output+"_bit",compare_z_score_bit,compare_norm_score_bit)
    writing(output,compare_z_score,compare_norm_score)
    