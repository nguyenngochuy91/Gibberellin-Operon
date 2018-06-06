#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : create operon tree
    Start   : 06/04/2017
    End     : 06/08/2017
'''
import os
import argparse
import shutil
from Bio import SeqIO
from ete3 import Tree
def parser_code():

    parser = argparse.ArgumentParser(description='The purpose of this script is to build a newick format gene tree from a list of genomes.')
    
    parser.add_argument("-G", "--genbank_directory", dest="genbank_directory", metavar="DIRECTORY", default='./db/',
                help="Folder containing all genbank files for use by the program.")
    
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="DIRECTORY", default='./operon_tree/',
                help="Directory where the results of this program will be stored.")
                
                
    
                                        
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
                target = info[1]
                max_percent = float(info[2])
                dic[f] = target               
                for line in infile.readlines():
                    info = line.split('\t')
                    target = info[1]
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
    print ("cmd1",cmd1)
    print ('Finish blasting:',target)

# use dir from best_gene , make fasta file
def make_target_fasta_operon(result, operon_dic,marker_fasta,marker_genes,has_blocks):
#    print operon_dic
#    print result
    out_fasta = open(marker_fasta,"w")
    out_fasta_dic = {} # key by the species, value will be combination of translation of each genes in the marker_genes
    for marker_gene in marker_genes:
        species_gene_dic = operon_dic[marker_gene]
        
        for file in result:
            name   = file.split('/')[-1].split(".ffc")[0]
            locus  = species_gene_dic[name]
            accession = locus.split('|')[0]
            if accession not in has_blocks:
                continue
#            print locus
#            print file
            for record in SeqIO.parse(file,"fasta"):
                if record.id == locus:
                    if accession not in out_fasta_dic:
                        out_fasta_dic[accession] = str(record.seq)
                    else:
                        out_fasta_dic[accession] += str(record.seq)
#    print out_fasta_dic
    for species in out_fasta_dic:
        out_fasta.write(">"+species+"\n" + out_fasta_dic[species] +"\n")
    out_fasta.close()
    
def make_newick_tree(marker_fasta, tree_outfile):
    ## using muscle
    # make the alignment file
    temp_align = '.'.join(marker_fasta.split('.')[:-1]) + '.aln' 
    cm1 ="muscle -in "+marker_fasta+ " -out "+temp_align +" -gapopen 0 "
    os.system(cm1)
    #make the tree using clustal
    cm2 ="clustalw -infile="+temp_align+" -tree=1"
    # have to wait for few second for the aln file actually comes out lol
    os.system(cm2)
    temp_tree = '.'.join(marker_fasta.split('.')[:-1]) + '.ph' # that's what this file gets named by default, and i'm sick of looking for the cmd line arg to fix.
    print (temp_tree)
    print ("modifying")
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
    has_blocks =[]
    infile = open('has_blocks.txt','r')
    for line in infile.readlines():
        has_blocks.append(line.strip())
    infile.close()
    parsed_args       = parser_code()
    genbank_directory = parsed_args.genbank_directory # db 

    output            = parsed_args.outfolder
    try:
        os.mkdir(output)
    except:
        print (output+" was already created")
    marker_genes       = ["KS","CPS","CYP117","SDR","FD","CYP114","CYP112"]
    # dic that store info for core genes in operon, key is the gene names, value is dictionary that key by species names, value 
    # are the locus info
    operon_dic = {}
    for marker_gene in marker_genes:
        try:
            os.mkdir(output+marker_gene)
        except:
            print (output+marker_gene +" was already created")
        result = [i for i in return_recursive_dir_files(genbank_directory) if i.split('/')[-1].split('.')[-1] == 'ffc']
        for file in result:
            name = file.split('/')[-1]
            parse_written(name,marker_gene+'.fa',output+marker_gene)
        # for each result of blast for a specie, get the locus of the best hit for gene
        # create a dic that map a genome to its best hit of gene
        operon_dic[marker_gene] = best_gene(output+marker_gene+'/')
    # write out a file .fa that store the translation gene for each species
    tree_dir          = output+'tree/'
    marker_fasta      = tree_dir + "distmat_marker.fa"
    tree_outfile      = tree_dir + "out_tree.nwk"
    try:
        os.mkdir(tree_dir)
    except:
        print (tree_dir+" was already created")
    make_target_fasta_operon(result, operon_dic, marker_fasta,marker_genes,has_blocks)
    # using the marker_fasta, create newick tree
    make_newick_tree(marker_fasta, tree_outfile)

main()