#!/usr/bin/python

# TODO: remove the last of the operon references.  this has not been done yet due to a piece of untested code. 
# until i get it teched out, i will not mess with teh variable names.
# This code needs to be cleaned up a lot... there is a lot of testing garbage that is no longer needed!
# It however, does work for the purposes of the paper submission.

from multiprocessing import Pool
import time
import os
import sys
#import simplejson as json
import json
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Copyright(C) 2015 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description='Convert a gene block file into a BLAST query using one or more reference organisms.')
    
    parser.add_argument("-i", "--infolder", dest="infolder", metavar="DIRECTORY", default='./genomes/',
                help="Folder containing all genbank files for use by the program. The refrence genbank file must reside here.")

                             
    parser.add_argument("-o", "--outfile", dest="outfile", metavar="FILE", default='./gene_block_query.fa',
                help="Resulting BLAST query file the is derived from a gene block file and refrence organism(s).")
                
    parser.add_argument("-b", "--gene_block_file", dest="gene_block_file", metavar="FILE", default='./regulonDB/gene_block_names_and_genes.txt',
                help="File which contains gene block information for use in gene block queries. The file format is gene_block_name followed by the constituent gene names, tab delineated.")
                
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")

    parser.add_argument("-r", "--refrence", dest="refrence", metavar="STRING", default = 'reference',
                help="Accession number of the refrence organism. This information is used to determine the product type of each gene (RNA/Protein), a necessary piece of information to classify the gene blocks that are under investigation.")
               
    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False,
                help="Suppresses most program text outputs.")
                  
    return parser.parse_args()
    
    
def check_options(parsed_args):   
    # I'm not issuing a message that hey, this file is there and will be overwritten, the program will just overwrite.
    # TODO: check that the folder that this is saved in exists
    
    if os.path.isdir(parsed_args.infolder):
        infolder = parsed_args.infolder
    else:
        print "The infolder %s does not exist." % parsed_args.infolder
        sys.exit()
    
    outfile = parsed_args.outfile
    
    # Check the gene block file
    if os.path.exists(parsed_args.gene_block_file):
        gene_block_file = parsed_args.gene_block_file
    else:
        print "The gene_block file %s does not exist." % parsed_args.gene_block_file
        sys.exit()
    
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)
    
    # I am not checking this option here, just leave it as default for the time being.
    refrence_list = [parsed_args.refrence]
    
    quiet = parsed_args.quiet
    
    return infolder, outfile, gene_block_file, num_proc, refrence_list, quiet


#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# convert a file of format : gene_block name then a list of the full names of the genes within that gene block
# into a dictionary that can be easily accessed or filtered later.
def parse_gene_block_file(fname):
    result = {}
    
    for line in [i.strip().split('\t') for i in open(fname).readlines()]:
        result.update({line[0]: line[1:]})
        
    return result

# This function creates a dictionary indexed by locus from the input genbank file
# and for my purposes now, it will index genes based on their annotation in genbank
# seq_type will allow us to determine the type of sequence returned in for the dict. default
# will be amino acid because this is a lot less noisy.
def return_genbank_dict(gb_file, key = 'locus', seq_type = 'amino_acid'):
    """Overview: This function will return a dictionary generated from a genbank file with key value supplied by caller.
       Returns: A dictionary created by the supplied genbank file (gb_file) indexed off the key value supplied.
       Default: The deafult key is locus, and this is generally the most useful key type since it is garanteed to be 
       unique within the genbank file. This condition is not necessarily true for any other attribute.
   """
    result = {}
    seq_record = SeqIO.parse(open(gb_file), "genbank").next()
    accession = seq_record.annotations['accessions'][0].split('.')[0]
    common_name = seq_record.annotations['organism'].replace(' ', '_')
    result.update({'accession': accession})
    result.update({'common_name': common_name})
    cnt = 0
    # loop over the genbank file
    unk_cnt = 1
    for fnum, feature in enumerate(seq_record.features):
        # here i simply check the gene coding type, and identify them in a way that can be used later.
        if feature.type == 'CDS' or feature.type == 'ncRNA' or feature.type == 'tRNA' or feature.type == 'mRNA' or feature.type == 'rRNA':
            start = feature.location.start
	       
            stop = feature.location.end
            #print start, stop
            strand = feature.strand
            synonyms = 'NONE'
            '''
            try: 
                gene = feature.qualifiers['gene'][0]
            except:
                gene = 'unknown'
            '''
            
            # this line might be wrong, just trying to get rid of an unneccessary try/except clause    
            if 'gene' in feature.qualifiers: 
                gene = feature.qualifiers['gene'][0]
            else:
                gene = 'unknown'
                    
            if 'gene_synonym' in feature.qualifiers:
                synonym_list = feature.qualifiers['gene_synonym'][0].replace(' ', '').split(';')
                synonyms = ':'.join(synonym_list)
            try:
                locus = feature.qualifiers['locus_tag'][0]
            except:
                try:
                    locus = feature.qualifiers['gene'][0]
                except:
                    locus = ''
                    print 'No locus associated. This should never be invoked meaning you are proper fracked. (The gbk file has an error).'
            try:
                seq = feature.qualifiers['translation']
                seq_type = 'Protein'
            except:
                cnt = cnt + 1
                seq = seq_record.seq[start:stop]
                seq_type = feature.type
                if feature.type == 'CDS':
                    seq_type = 'Pseudo_Gene'
                    # attempt to fix an error
                    if strand == 1:
                        seq = seq.translate()
                    else:
                        seq = seq.reverse_complement().translate()
            gc = "%2.1f" % GC(seq_record.seq[start:stop])
            # Debugging something odd
            
                #print feature.qualifiers['gene_synonym']
            #method = "exact"
            if key == 'locus':
                if gene == 'unknown':
                    new_gene = 'unknown_' + str(unk_cnt)
                    header = '|'.join([accession, common_name, locus, new_gene, str(start), str(stop), str(strand), seq_type, synonyms, gc])
                    result.update({locus: [header, ''.join(seq)]})
                    unk_cnt +=1
                else:
                    header = '|'.join([accession, common_name, locus, gene, str(start), str(stop), str(strand), seq_type, synonyms, gc])
                    result.update({locus: [header, ''.join(seq)]})
                    try:
                        for syn in synonym_list:
                            result.update({syn: [header, ''.join(seq)]})
                    except:
                        pass

                # result.update({locus: (locus, gene, seq, seq_type, synonyms)})
            elif key == 'annotation':
                if gene == 'unknown':
                    new_gene = 'unknown_' + str(unk_cnt)
                    header = '|'.join([accession, common_name, locus, new_gene, str(start), str(stop), str(strand), seq_type, synonyms, gc])
                    result.update({new_gene: [header, ''.join(seq)]})
                    unk_cnt +=1
                else:
                    header = '|'.join([accession, common_name, locus, gene, str(start), str(stop), str(strand), seq_type, synonyms, gc])
                    result.update({gene: [header, ''.join(seq)]})
                    try:
                        for syn in synonym_list:
                            result.update({syn: [header, ''.join(seq)]})
                    except:
                        pass

    #print 'The number of non-protein regions in %s is: %i.' % (common_name, cnt)
    return result


# This function will allow me to do the main work of make_gene_block_fasta, but allow parallel
# processing. I will have to make a parallel array, which will take some time to learn.  i 
# will keep this stub for later to implement. 
def parallel_gene_block_fasta(genome):
    organism = genome.split('/')[-1].split('.')[0]
    organism_dict_for_recovery = {}
    org_dict = return_genbank_dict(genome)
    organism_dict_for_recovery.update({organism: org_dict})
    return (organism, org_dict)



# This function will make a BLAST query from the parsed gene_block file.  The default behavior of the function is 
# to make a query file from all genes in all organisms that are annotated. Later the results will be sorted based 
# the needs of the programmer. The defaulted variables allow a single gene_block to be chosen individually.  The program
# will also store the results of this function in a folder titled blast_query_files, in the recovery folder.
def make_gene_block_fasta2(gene_list, genbank_list, num_processors, outfile, ref_list):
    
    refrence_file_path_list = [i for i in genbank_list if os.path.basename(i).split('.')[0] in ref_list]
    #print "refrence_file_path_list", refrence_file_path_list
    pool = Pool(processes = num_processors)
    organism_dict_for_recovery = dict(pool.map(parallel_gene_block_fasta, refrence_file_path_list))
    
    protein_match = []
    rna_match = []
    pseudogene_match = []
    missing_list = []
    
    refrence_prot = []
    refrence_rna = []

    # This list should be updated should other types of RNAs be annotated later as parts of (siRNA comes to mind).
    RNA_codes = ['rRNA', 'tRNA', 'ncRNA']
    for org in organism_dict_for_recovery.keys():
        for gene in gene_list:
            if gene in organism_dict_for_recovery[org].keys():
                if organism_dict_for_recovery[org][gene][0].split('|')[7] == 'Protein':
                    outseq = SeqRecord(Seq(organism_dict_for_recovery[org][gene][1]),
                    id=organism_dict_for_recovery[org][gene][0], description = '')
                    protein_match.append(outseq)
                    if org in ref_list:
                        outseq = SeqRecord(Seq(organism_dict_for_recovery[org][gene][1]),
                        id=organism_dict_for_recovery[org][gene][0], description = '')
                        refrence_prot.append(outseq)
                elif organism_dict_for_recovery[org][gene][0].split('|')[7] in RNA_codes:
                    outseq = SeqRecord(Seq(organism_dict_for_recovery[org][gene][1]),
                    id=organism_dict_for_recovery[org][gene][0], description = '')
                    rna_match.append(outseq)
                    if org in ref_list:
                        outseq = SeqRecord(Seq(organism_dict_for_recovery[org][gene][1]),
                        id=organism_dict_for_recovery[org][gene][0], description = '')
                        refrence_rna.append(outseq)
                elif organism_dict_for_recovery[org][gene][0].split('|')[7] == 'Pseudo_Gene':
                    outseq = SeqRecord(Seq(organism_dict_for_recovery[org][gene][1]),
                    id=organism_dict_for_recovery[org][gene][0], description = '')
                    pseudogene_match.append(outseq)
                    if org in ref_list:
                        outseq = SeqRecord(Seq(organism_dict_for_recovery[org][gene][1]),
                        id=organism_dict_for_recovery[org][gene][0], description = '')
                        refrence_prot.append(outseq)
                else:
                    print "There was an error in function make_gene_block_fasta2 from script make_operon_query.py in",  organism_dict_for_recovery[org][gene][0]
            else: # The gene is missing, we will look for it at a later stage in the program
                item = '\t'.join([org, gene])
                missing_list.append(item)
            

        
    handle = open(outfile, 'w')       
    SeqIO.write(protein_match, handle,"fasta")
    handle.close()
    
    '''    
    handle = open(folder + 'rna_matches.fa', 'w')       
    SeqIO.write(rna_match, handle,"fasta")
    handle.close()
        
    handle = open(folder + 'pseudogene_matches.fa', 'w')       
    SeqIO.write(pseudogene_match, handle,"fasta")
    handle.close()
        
    handle = open(folder + 'missing_gene_block_genes.txt', 'w')
    handle.write('\n'.join(missing_list))
    handle.close()
    
    
    # The goal here is to return the refrence fasta files, so that we can do self homolog dict, which is needed to detect
    # fusion proteins later.
    #ref_prot_outfile = folder + 'refrence_prot.fa'
    #ref_rna_outfile = folder + 'refrence_rna.fa'
    
        
    #print "len(gene_list)", len(gene_list)
    handle = open(ref_prot_outfile, 'w')
    SeqIO.write(refrence_prot, handle, "fasta")
    #print "len(refrence_prot)", len(refrence_prot)
    handle.close()
        
    handle = open(folder + 'refrence_rna.fa', 'w')
    SeqIO.write(refrence_rna, ref_rna_outfile, "fasta")
    #print "len(refrence_rna)", len(refrence_rna)
    handle.close()
    

    return ref_prot_outfile, ref_rna_outfile, organism_dict_for_recovery
    '''

# This function validates the gene_blocks in the references. any gene_blocks which cannot be located in their entirety are omitted
# it returns two dicts.  one has just the gene_blocks and a list of their genes, the other has gene_blocks, list of genes with the
# type of product produced
def categorize_gene_blocks(ref_org, genbank_list, gene_block_dict):
    ref_path = [i for i in genbank_list if i.split('/')[-1].split('.')[0] in ref_org]
    print "ref_org",ref_org
    print "ref_path",ref_path
    # dict of the gene_blocks that are validated (could cind all the constituent genes in the reference)
    result = {}
    # dict of the gene_blocks that are validated (could cind all the constituent genes in the reference) and the type of product
    result_categorized = {}
    
    #print "gene_block_dict", len(gene_block_dict), gene_block_dict
    
    #ref_dict = return_genbank_dict(ref_path)
    # basically need to update this for a list of refrence organisms... this is not done just yet.
    ref_dict = return_genbank_dict(ref_path[0])
     
    for gene_block in sorted(gene_block_dict.keys()):
        gene_block_error = False
        type_list = []
        print "gene_block_dict",gene_block_dict[gene_block]
        for gene in gene_block_dict[gene_block]:
            print "gene",gene
            try:
                gene_type = ref_dict[gene][0].split('|')[7]
                # if protein
                if gene_type == 'Protein' or gene_type == 'Pseudogene':
                    type_list.append("%s:p" % gene)
                # if RNA
                else:
                    type_list.append("%s:r" % gene)
                    #print "RNA", gene, gene_type

            except:
                print "gene_block", gene_block, "is not usable as an gene_block, I should remove it. The gene in error is in", gene
                gene_block_error = True

        # no errors
        if not gene_block_error: 
            result.update({gene_block: gene_block_dict[gene_block]})
            result_categorized.update({gene_block:type_list})
            #print "result", result
            #print "result_categorized", result_categorized

    return result, result_categorized




def main():
    
    start = time.time()

    parsed_args = parser_code()
    
    infolder, outfile, gene_block_file, num_proc, refrence_list, quiet = check_options(parsed_args)
    
    if not quiet:
        print infolder, outfile, gene_block_file, num_proc, refrence_list, quiet
    
    # This section of code will return a parsed gene block file as a dictionary keyed by gene block name
    #parsed_gene_block_file = './regulonDB/gene_block_names_and_genes_unfiltered.txt'
    # later, this should be done through the cmd line interface... just running out of time here

    gene_block_dict = parse_gene_block_file(gene_block_file)

    # This section of code will return the full pathways to the genbank files of interest
    genbank_list = returnRecursiveDirFiles(infolder)
    
    # check to make sure that all gene blocks that we are looking at can be found in their entirety before we use them as 
    # part of the data set.  The validated gene block dict contains information about the type (protein/RNA) of gene product.
    validated_gene_block_dict, validated_gene_block_dict_more_info = categorize_gene_blocks(refrence_list, genbank_list, gene_block_dict)

    # this is a list of all genes that we have in the dataset, reguardless of wether we can find them in the reference or not.
    unvalidated_gene_list = []
    
    # this is a list of all the genes that are in the validated dataset.  we can find all the genes in the reference organism for the entire operon.
    validated_gene_list = []

    for gene_block in gene_block_dict.keys():
        unvalidated_gene_list = unvalidated_gene_list + gene_block_dict[gene_block]
        
    for gene_block in validated_gene_block_dict.keys():
        validated_gene_list = validated_gene_list + validated_gene_block_dict[gene_block]

    #ref_prot_outfile, ref_rna_outfile, org_annotation_dict = make_gene_block_fasta2(validated_gene_list, genbank_list, num_proc, outfolder, refrence)
    make_gene_block_fasta2(validated_gene_list, genbank_list, num_proc, outfile, refrence_list)
    
    # currently, we are only looking at protein sequences. The following code reflects this.
    '''
    self_homolog_dict = return_self_homolog_dict(validated_operon_dict, ref_prot_outfile, outfolder + 'operon_homolog_dict.json', outfolder)
    
    # at this point i have the self homolog dict, i have a file that has all the e.coli operons validated (for one version of K-12)
    # and i have the operons validated per operon, at least in terms of prot or RNA.  now what i am going to do it to take the operons 
    # individually, make a folder for them in the outfolder, named for the operon, and create a fasta file of all the annotated examples.
    # once this is done, then we will run CD hit on the individual genes to get a set of represenatives per gene.  the clustering will be
    # done on the a.a. seq at 70%.  Then i would like to screne the clusters to try to remove possible mis-annotations.
    
    make_gene_block_individual_gene_fasta(validated_operon_dict, org_annotation_dict, refrence, outfolder)

    '''
    # easy way to run this, using all the defaults that make sense
    # ./make_operon_query.py -i /home/dave/Desktop/all_genbank -o ./operon_query.fa -p ./regulonDB/operon_names_and_genes.txt
    
    # ./make_operon_query.py -i ./genomes/ 

    if not quiet:
        print time.time() - start
    
if __name__ == '__main__':
    main()
