#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : boostrapping
    Start   : 06/04/2017
    End     : 06/08/2017
'''

from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.Consensus import *
from Bio import AlignIO,Phylo
import argparse
import time
from ete3 import *
def parser_code():

    parser = argparse.ArgumentParser(description='The purpose of this script is generate branch support of our generated tree.')
    
    parser.add_argument("-t", "--input_tree", help="Input newick tree")
    
    parser.add_argument("-i", "--alignment", help="Multiple alighment file")
    
    parser.add_argument("-n","--number",help="How many tree to draw")
    
    parser.add_argument("-o", "--ouput_tree",
                help="Output support tree")
    parser.add_argument("-f", "--filter", dest="filter", metavar="FILE", default='NONE',
                help="File restrictiong which accession numbers this script will process. If no file is provided, filtering is not performed.")
    
                                        
    return parser.parse_args()
def filterMSA(msa,filter_list,):
    output = open("filter_msa.fasta","w")
    msa          = AlignIO.read(msa,"fasta")
    for item in msa:
        if item.id in filter_list or item.name in filter_list:
            output.write(">{}\n{}\n".format(item.id,item.seq))
    output.close()
    
if __name__ == '__main__':
    arguments    = parser_code()
    filter_file  = arguments.filter
    number       = int(arguments.number)
    input_tree   = Phylo.read(arguments.input_tree,"newick")
    msa          = arguments.alignment
    output_tree  = arguments.ouput_tree
    if filter_file!='NONE':
        filter_list = [i.strip() for i in open(filter_file,"r")]
        filterMSA(msa,filter_list)
        msa = AlignIO.read("filter_msa.fasta","fasta")
    else:
        msa          = AlignIO.read(msa,"fasta")
    calculator   = DistanceCalculator('blosum62')
    constructor  = DistanceTreeConstructor(calculator)
    start        = time.time()
    trees        = bootstrap_trees(msa, number, constructor)
    trees        = list(trees)
    # set "KE136308.1" as outgroup
    for t in trees:
        t.root_with_outgroup("KE136308.1")
    support_tree = get_support(input_tree, trees)
    print (support_tree)
    Phylo.write(support_tree,output_tree+".nwk","newick")
#    tree         = Tree(output_tree+".nwk")
#    tree.render(output_tree+".pdf")
    stop  = time.time()
    print (stop-start)
    tree = Tree(output_tree+".nwk")
    tree_style = TreeStyle()
    tree_style.min_leaf_separation = 5
    tree_style.extra_branch_line_type = 0
    tree_style.draw_guiding_lines=True
    tree_style.guiding_lines_type = 1
    for node in tree.iter_descendants("postorder"):
        if node.is_leaf():
            node.add_face(TextFace(node.name,fgcolor = 'blue'), column =0, position ="aligned")
        else:
            node.add_face(TextFace(node.support),0,"branch-top")
    tree.render(output_tree+'.pdf',dpi=1000,tree_style=tree_style)