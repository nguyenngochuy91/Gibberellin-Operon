#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : blasting each 
    Start   : 09/28/2016
    End     : /2016
'''
import os
import argparse
import shutil
from ete3 import Tree
from Bio import SeqIO
def parser_code():

    parser = argparse.ArgumentParser(description='The purpose of this script is to build a newick format gene tree from a list of genomes.')
    
    parser.add_argument("-i", "--genbank_directory", dest="genbank_directory", metavar="DIRECTORY", default='./db/',
                help="Folder containing all genbank files for use by the program.")
    
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="DIRECTORY", default='./gene_tree/',
                help="Directory where the results of this program will be stored.")
                
                 
    parser.add_argument("-m", "--marker_gene", dest="marker_gene", metavar="STRING", default='rpob',
                help="This is a single marker gene that will be used to construct phylogenetic trees (CPS, CYP117,CYP114).")
                
    
                                        
    return parser.parse_args()
    
# get the files inside a directory
def return_recursive_dir_files(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            if '.DS_Store' not in f:
                fname = os.path.join(path, f)
                if os.path.isfile(fname):
                    result.append(fname)
    return result

# blast each full sequence against the traslation of gene in gilbberin

# create dic using all the blast file in the gene_tree/gene dir
def best_gene(root_dir):
    dic = {}
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            if '.DS_Store' not in f:
                fname = os.path.join(path, f)
                infile = open(fname,'r')
                line = infile.readline()
                info = line.split('\t')
                if len(info)<=1:
                    continue
                target = info[1]
                max_percent = float(info[2])
                dic[f] = target               
                for line in infile.readlines():
                    info = line.split('\t')
                    query = info[0].split("|")
                    gene_length      = (int(query[5])-int(query[4]))/3
                    alignmentLength  = float(line[3]) 
                    if alignmentLength/gene_length>=(2/3.0):   
                        temp = target+"|False"
                        target= temp
                    percent = float(info[2])
                    if percent> max_percent:
                    # get the one with highes percentage
                        dic[f] = target
                        max_percent = percent
    return dic
                            
# blast gene vs db database              
def parse_written(target,query,blast_dir):
    cmd1 = ('blastp -query '+query+' -outfmt 6 -out '+blast_dir+'/'+target.split('.ffc')[0]+
        ' -subject db/'+target+' -num_threads 8')
    os.system(cmd1)
    #print ("cmd1",cmd1)
    print ('Finish blasting:',target)

# use dir from best_gene , make fasta file
def make_target_fasta(result, species_gene_dic,marker_fasta,has_blocks,has_CYP114):
    out_fasta = open(marker_fasta,"w")
    out_fasta_dic = {}
    for file in result:
        name   = file.split('/')[-1].split(".ffc")[0]
        locus  = species_gene_dic[name]
        infile = open(file,"r")
        line   = infile.readline()
        info = line.split('|')
        if info[0][1:] in has_blocks and info[0][1:] in has_CYP114:
            accession = locus.split('|')[0]
            for record in SeqIO.parse(file,"fasta"):
                if record.id == locus:
                    out_fasta_dic[accession] = str(record.seq)
    for species in out_fasta_dic: 
        out_fasta.write(">"+species+"\n" + out_fasta_dic[species] +"\n")
    out_fasta.close()
    
def make_newick_tree(marker_fasta, tree_outfile):
    ## using muscle
    # make the alignment file
    temp_align = '.'.join(marker_fasta.split('.')[:-1]) + '.aln' 
    cm1 ="muscle -in "+marker_fasta+ " -out "+temp_align
    os.system(cm1)
    #make the tree using clustal
    cm2 ="clustalw -infile="+temp_align+" -tree=1"
    # have to wait for few second for the aln file actually comes out lol
    os.system(cm2)
    temp_tree = '.'.join(marker_fasta.split('.')[:-1]) + '.ph' # that's what this file gets named by default, and i'm sick of looking for the cmd line arg to fix.
    #modify for negative branch
    modify_tree = '.'.join(marker_fasta.split('.')[:-1]) + '.new'
    cm3 = "sed -e 's,:-[0-9\.]\+,:0.0,g' "+temp_tree+" > "+modify_tree   
    os.system(cm3)
    tree = Tree(modify_tree)
    tree.set_outgroup("KE136308.1")
    # dealing with negative branch length
    #print "marker_fasta",marker_fasta
    #print "temp_tree", temp_tree
    # move the created tree file to the location i say its going
    tree.write(outfile= tree_outfile)
    
def main():
    try:
        os.mkdir("gene_tree")
    except:
        print ("gene tree dir is already created")
    has_blocks =[]
    infile = open('has_block.txt','r')
    for line in infile.readlines():
        has_blocks.append(line.strip())
        
    has_CYP114= []
    CYP114     = open("CYP114.txt","r")
    for line in CYP114.readlines():
        has_CYP114.append(line.strip())    
    infile.close()
    parsed_args       = parser_code()
    genbank_directory = parsed_args.genbank_directory # db 
    marker_gene       = parsed_args.marker_gene # get the gene name
    output            = parsed_args.outfolder
    tree_dir          = output+marker_gene+'/tree/'
    marker_fasta      = tree_dir + "distmat_marker.fa"
    tree_outfile      = tree_dir + "out_tree.nwk"
    try:
        os.mkdir(output+marker_gene)
    except:
        print (output+marker_gene+" already created")
    result = [i for i in return_recursive_dir_files(genbank_directory) if i.split('/')[-1].split('.')[-1] == 'ffc']
    for file in result:
        name = file.split('/')[-1]
        parse_written(name,marker_gene+'.fa',output+marker_gene)
    # for each result of blast for a specie, get the locus of the best hit for gene
    # create a dic that map a genome to its best hit of gene
    species_gene_dic = best_gene(output+marker_gene+'/')
#    print species_gene_dic
    # write out a file .fa that store the translation gene for each species
    try:
        os.mkdir(tree_dir)
    except:
        print (tree_dir+" already created")
    make_target_fasta(result, species_gene_dic, marker_fasta,has_blocks,has_CYP114)
    # using the marker_fasta, create newick tree
    make_newick_tree(marker_fasta, tree_outfile)

main()