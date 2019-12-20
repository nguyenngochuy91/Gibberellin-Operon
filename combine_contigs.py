#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : for each file in family genome folder, string all contigs into
             1 big pseudo genome
    Start   : 09/28/2016
    End     : /2016
'''
import os
import time
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
from Bio.SeqRecord import SeqRecord
# arguments
def get_arguments():
    parser = argparse.ArgumentParser(description='combine all contigs into 1 file, and provide file to describe.')
    parser.add_argument("-G", "--fasta_directory", default='./family_contigs_folder/',
                help="Folder containing all fasta files for use by the program.(/family_contigs_folder/)")
    parser.add_argument("-o", "--outfolder_directory", default='./family_genome_folder/',
                help="Directory where the combination of contigs will be stored(/family_genome_folder/)")
    parser.add_argument("-m", "--contig_map_directory", default='./contigs_map/',
                help="File that store the index of each contigs of the contigs combining process(/contigs_map/)")
    parser.add_argument("-t", "--threshold", default=10000,
                help="Threshold for each contig to be kept")                
    return parser.parse_args()

#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def return_recursive_dir_files(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            if '.DS_Store' not in f:
                fname = os.path.join(path, f)
                if os.path.isfile(fname):
                    result.append(fname)
    return result
'''
    function: giving a fasta file of contigs, combine the sequence and output into
    a fasta format using the description of the first contig. Also a dictionary that map
    the index of each contig in the combination. In addition, also filter those that don't have
    good quality
    input : fasta (in_fasta), 
    output: dic (out_fasta_dic), dic (out_map_dic)
'''
def combine_sequence(in_fasta,threshold):
    record_dict = SeqIO.index(in_fasta, "fasta") # index the record
    # initiate an empty Sequence string
    combined_string = Seq("",SingleLetterAlphabet())
    # fasta_dic
    out_fasta_dic ={} 
    # dictionary map index of each contig
    out_map_dic ={}
    start = 0
    count = 0
    for record in sorted(record_dict):
        count +=1 
        new_contig = record_dict[record].seq
        if len(new_contig) >= threshold:
            length = len(new_contig)
            combined_string += new_contig
            out_map_dic[record_dict[record].id] = [start,start+length]
            start += length # increment the start position
        if count == 1:
            out_fasta_dic["id"] = record_dict[record].id
            out_fasta_dic["description"] = record_dict[record].description
    out_fasta_dic["sequence"] = combined_string
    return out_fasta_dic,out_map_dic
'''
    function: giving out_fasta, and out_map, write this into a file
    input : dic (out_fasta_dic), dic (out_map)
    output: fasta (out_fasta), text (out_map)
'''
def write_file(out_fasta_dic,out_map_dic,out_fasta,out_map):
    # write the record from out_fasta_dic to out_fasta 
    record = SeqRecord(out_fasta_dic["sequence"], id = out_fasta_dic["id"],
                   description=out_fasta_dic["description"])
    SeqIO.write(record, out_fasta, "fasta") 
    # write the out_map
    outfile = open(out_map,'w')
    for contigs in out_map_dic:
        info  = out_map_dic[contigs]
        start = info[0]
        stop  = info[1]
        outfile.write(contigs +':'+str(start)+','+str(stop)+'\n')
    outfile.close()

    
if __name__ == "__main__":
    args                    = get_arguments()
    fasta_directory         = args.fasta_directory
    outfolder_directory     = args.outfolder_directory
    contig_map_directory    = args.contig_map_directory
    threshold               = int(args.threshold)
    
    # make the dir
    if not os.path.isdir(outfolder_directory):
        os.mkdir(outfolder_directory)
    if not os.path.isdir(contig_map_directory):
        os.mkdir(contig_map_directory)
    # walk the fasta_directory
    start  = time.time()
    result = return_recursive_dir_files(fasta_directory)
    for path in result:
        name = path.split('/')[-1]
        out_fasta_dic,out_map_dic = combine_sequence(path,threshold)
        if len(out_fasta_dic["sequence"]) > 0:
            out_fasta = outfolder_directory + name 
            out_map   = contig_map_directory + name
            write_file(out_fasta_dic,out_map_dic,out_fasta,out_map)
            print ("Done with this file:",name)
    stop = time.time()
    print ("Eh, it took a while",stop-start)