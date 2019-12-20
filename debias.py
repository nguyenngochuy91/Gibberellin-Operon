#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : debias
    Start   : 06/04/2017
    End     : 06/08/2017
'''
import os
import argparse
import shutil
from Bio import SeqIO
from ete3 import Tree
def parser_code():

    parser = argparse.ArgumentParser(description='The purpose of this script to debias tree based on parameter')
    
    parser.add_argument("-i", "--input_tree", help="Input tree that we want to debias")
    
    parser.add_argument("-o", "--output_tree", help="Output tree to be store.")
                
    parser.add_argument("-s", "--tree_size", help="Size of the tree for output, for example reduced from 149 leaves to 100 leaves")
 
    parser.add_argument("-k", "--keeper", help="Force to include the following speciese")     
    parser.add_argument("-a", "--accession", help="Accession file accession_to_common.txt")                                            
    return parser.parse_args()
    
def parse_pda(handle):
    for line in handle.readlines():
        if len(line) > 100:
            return line.strip()
    
    
if __name__ == "__main__":
    args  = parser_code()
    input_tree = args.input_tree
    output_tree = args.output_tree
    accession_file = args.accession
    difference = '/'.join(input_tree.split('/')[:-1])+"/not_included.txt"
    filter_file = '/'.join(input_tree.split('/')[:-1])+"/filter.txt"
    t = Tree(input_tree, format = 1)
    leaves_full = t.get_leaf_names()
    size = int(args.tree_size)
    keep = args.keeper
    cmd1 = "./pda -k {} {} {} -if {}".format(size,input_tree,output_tree,keep)
    os.system(cmd1)
    tree = parse_pda(open(output_tree,"r"))
    t = Tree(tree,format = 1)
    leaves_partial = t.get_leaf_names()
    print (leaves_partial)
    not_included = [item for item in leaves_full if item not in leaves_partial]
    t.set_outgroup("KE136308.1")
    t.write(outfile = output_tree, format = 1)
    outfile = open(difference,'w')
    for item in not_included:
        outfile.write(item+"\n")
    outfile.close()
    outfile = open(filter_file,"w")
    infile = open(accession_file,"r")
    for line in infile.readlines():
        for leaf in leaves_partial:
            if leaf in line:
                outfile.write(line)
                break
    outfile.close()
    infile.close()
        