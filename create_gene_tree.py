#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : blasting each 
    Start   : 09/28/2016
    End     : /2016
'''
import os
import argparse
import shutil
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
                for line in infile.readlines():
                    line = line.split('\t')[1]
                    line = line.split('|')[2]
                    dic[f] = line
    return dic
                            
# blast gene vs db database              
def parse_written(target,query,blast_dir):
    cmd1 = ('blastp -query '+query+' -outfmt 6 -out '+blast_dir+'/'+target.split('.ffc')[0]+
        ' -subject db/'+target+' -num_threads 8 -max_target_seqs 1')
    os.system(cmd1)
    #print ("cmd1",cmd1)
    print ('Finish blasting:',target)

# use dir from best_gene , make fasta file
def make_target_fasta(result, species_gene_dic,marker_fasta,has_blocks):
    out_fasta = open(marker_fasta,"w")
    for file in result:
        name   = file.split('/')[-1].split(".ffc")[0]
        locus  = species_gene_dic[name]
        infile = open(file,"r")
        line   = infile.readline()

        while len(line) > 0:
            if line[0] ==">": 
                count = 0
                info = line.split('|')
                if info[0][1:] not in has_blocks:
                    print (info[0][1:])
                    break
                accesion = info[0]+"\n"
                current_locus = info[2]
                if locus == current_locus : # write this into out_fasta
                    out_fasta.write(accesion) # write the header
                    string =""
                    next_line = infile.readline().strip()
                    while next_line[0] != ">":
                        string += next_line
                        next_line = infile.readline().strip()
                        count +=1
                        if count >30:
                            break
                    string += "\n"
                    out_fasta.write(string)
                    break          # break when found all           
            line = infile.readline()
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
    print temp_tree
    print "modifying"
    #modify for negative branch
    modify_tree = '.'.join(marker_fasta.split('.')[:-1]) + '.new'
    cm3 = "sed -e 's,:-[0-9\.]\+,:0.0,g' "+temp_tree+" > "+modify_tree   
    os.system(cm3)
    # dealing with negative branch length
    #print "marker_fasta",marker_fasta
    #print "temp_tree", temp_tree
    # move the created tree file to the location i say its going
    shutil.copy(modify_tree, tree_outfile)
    
def main():
    try:
        os.mkdir("gene_tree")
    except:
        print ("gene tree dir is already created")
    has_blocks =[]
    infile = open('has_blocks.txt','r')
    for line in infile.readlines():
        has_blocks.append(line.strip())
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
        print output+marker_gene+" already created"
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
        print tree_dir+" already created"
    make_target_fasta(result, species_gene_dic, marker_fasta,has_blocks)
    # using the marker_fasta, create newick tree
    make_newick_tree(marker_fasta, tree_outfile)

main()