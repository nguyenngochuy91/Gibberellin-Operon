#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : create operon tree  using distance from operon
                given the file contains orthologs for each species, count the edit distance
                for each events for each pair of species. Then, normalize them and provide a distance
                between any 2 species that can be used to build operon trees.
    Start   : 06/04/2017
    End     : 06/08/2017
'''
import os
import argparse
import csv
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor,_DistanceMatrix
from Bio import Phylo
def parser_code():

    parser = argparse.ArgumentParser(description='The purpose of this script is to build a newick format gene tree from a list of genomes.')
    
    parser.add_argument("-i", "--input",default='ancestral_reconstruction/new_result/gibberellin',
                help="converted result")
    
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="DIRECTORY", default='./operon_tree_distance/',
                help="Directory where the results of this program will be stored.")
                
                
    
                                        
    return parser.parse_args()
###############################################################################
## helpers to calculate distance
    
# get the gene set from blocks
def get_gene_set(blocks):
    genes = set()
    for block in blocks:
        for gene in block:
            genes.add(gene)
    return genes
# get the duplication genes from blocks if available
def get_duplication_set(blocks):
    duplications = set()
    for block in blocks:
        dic = {}
        for gene in block:
            if gene not in dic:
                dic[gene]=1
            else:
                duplications.add(gene)
    return duplications

# provide the relevant blocks and calculate the difference in size of 2 given orthoblocks
def get_split_difference(blocks1,blocks2,intersection):
    new_blocks1 =[]
    new_blocks2 = []
    for item in blocks1:
        new_block = ""
        for gene in item:
            if gene in intersection:
                new_block+=gene
        if new_block!= "":
            new_blocks1.append(new_block)
    for item in blocks2:
        new_block = ""
        for gene in item:
            if gene in intersection:
                new_block+=gene
        if new_block!= "":
            new_blocks2.append(new_block)
    return abs(len(new_blocks2)-len(new_blocks1))


###############################################################################
## parsing file, creating matrix distance, creating trees
# get the orthologs
def parse_operon(infile,has_blocks):
    dic ={}
    infile.readline()
    for line in infile.readlines():
        line  = line.strip().split(':')
        species = line[0]
        if species not in has_blocks:
            continue
        ortholog = line[1]
        # removing p and k in ortholog
        ortholog = ortholog.replace('p','')
        ortholog =  ortholog.replace('k','')
        blocks = ortholog.split('|')
        ortholog = '|'.join([item for item in blocks if item!=''])
        genes    = get_gene_set(blocks)
        duplications = get_duplication_set(blocks)
        dic[species] = [ortholog,genes,duplications]
    return dic
    
# normalize a given list
def normalize(my_list,my_min,my_max):
    for item in my_list[1:]:
        for i in range(1,len(item)):
            item[i] = round((item[i]-my_min)/float(my_max-my_min),4)
            

# write from list into csv
def writing(filename,matrix1):
    writer=csv.writer(open(filename + ".csv","w"))
    for item in matrix1:
        writer.writerow(item)
    
"""
    function: given a dictionary, key is species, value is [ortholog, set of genes, set of duplications gebes]
                provide the 3 matrix of distance for each events for each pairs of species. In addition, it will be normalize
                and then add uptogether for the total sum, then normalize again to provide the final distance
    input   : dic
    output  :  list
"""
def get_total_distance(dic):
    names          = sorted(dic)
    deletion =[['']]
    dupication =[['']]
    split = [['']]
    total = [['']]
    
    deletion[0].extend(names)
    dupication[0].extend(names)
    split[0].extend(names)
    total[0].extend(names)
    # calculating all events distance for each pair of species
    for i in range(len(names)):
        deletion_current = [names[i]]
        dupication_current = [names[i]]
        split_current = [names[i]]
        for j in range(len(names)):
            # get the info from our dic given the species name:
            info_i = dic[names[i]]
            info_j = dic[names[j]]
            # if the 2 orthologs are the same, then just set all value to 0 for the 3 events
            if info_i[0] == info_j[0]:
                deletion_current.append(0)
                dupication_current.append(0)
                split_current.append(0)
                
                continue
            # if they are not equal, calculate the difference
            ## deletion event            
            deletion_difference = info_i[1].symmetric_difference(info_j[1])
            deletion_current.append(len(deletion_difference))
            ## duplication event
            dupication_difference = info_i[2].symmetric_difference(info_j[2])
            dupication_current.append(len(dupication_difference))           
            ## split event
            intersection  = info_i[1].intersection(info_i[2])
            split_difference = get_split_difference(info_i[0],info_j[0],intersection)
            split_current.append(split_difference)
        deletion.append(deletion_current)
        dupication.append(dupication_current)
        split.append(split_current)
    # write out csv file for debugging purpose
    writing("deletion_matrix",deletion)
    writing("dupication_matrix",dupication)
    writing("split_matrix",split)
    # store the max min info into dictionary, then normalize them

    deletion_max = max([max(item[1:]) for item in deletion[1:]])
    deletion_min = min([min(item[1:]) for item in deletion[1:]])
    normalize(deletion,deletion_min,deletion_max)
    
    dupication_max = max([max(item[1:]) for item in dupication[1:]])
    dupication_min = min([min(item[1:]) for item in dupication[1:]])
    normalize(dupication,dupication_min,dupication_max)
    
    split_max = max([max(item[1:]) for item in split[1:]])
    split_min = min([min(item[1:]) for item in split[1:]])
    normalize(split,split_min,split_max)
    # adding all these 3 events, and normalize them
    for i in range(len(names)):
        total_current =[names[i]]
        for j in range(len(names)):
            total_current.append(round(((deletion[i+1][j+1]+dupication[i+1][j+1]+split[i+1][j+1])),4))
        total.append(total_current)
    writing("total_matrix",total)   
    return total
"""
    function: given a distance matrix, construction the tree
    input   : Distance matrix from phylo
    output  : nwk tree 
"""
def create_matrix_distance(dic):
    names = sorted(dic)
    full_matrix = [item[1:] for item in get_total_distance(dic)[1:]]
    matrix = [full_matrix[i][:i+1] for i in range(len(names))]
    dm = _DistanceMatrix(names,matrix)
    return dm
def make_newick_tree(dm):
    constructor = DistanceTreeConstructor()
    upgmatree = constructor.upgma(dm)
    njtree = constructor.nj(dm)
    upgmatree.root_with_outgroup({'name': "KE136308.1"})
    njtree.root_with_outgroup({'name': "KE136308.1"})
    return upgmatree,njtree
def main():
    has_blocks =[]
    infile = open('has_blocks.txt','r')
    for line in infile.readlines():
        has_blocks.append(line.strip())
    infile.close()
    parsed_args       = parser_code()
    operon = parsed_args.input # matrix file
    outfolder       = parsed_args.outfolder
    dm = create_matrix_distance(parse_operon(open(operon,'r'),has_blocks))
#    print dm
    try:
        os.mkdir(outfolder)
    except:
        print outfolder +" already created"
    upgmatree,njtree = make_newick_tree(dm)
    Phylo.write(upgmatree,outfolder+"/upg_tree.nwk","newick")
    cm3 = "sed -e 's,:-[0-9\.]\+,:0.0,g' "+outfolder+"/upg_tree.nwk"+" > "+outfolder+"/new_upg.nwk" 
    os.system(cm3)
    Phylo.write(njtree,outfolder+"/nj_tree.nwk","newick")
    cm4 = "sed -e 's,:-[0-9\.]\+,:0.0,g' "+outfolder+"/nj_tree.nwk"+" > "+outfolder+"/new_nj.nwk" 
    os.system(cm4)
main()