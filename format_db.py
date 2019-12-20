#!/usr/bin/env python

from multiprocessing import Pool
import time
import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

# Copyright(C) 2015 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment


# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description='Convert all genbank files found in a specified folder, and optionally a file containing the accessions that you wish to include, and create BLAST searchable databases from them.')

    #parser.add_argument("-i", "--infolder", dest="infolder", metavar="DIRECTORY", default='./genomes/',
    #            help="Folder containing all genbank files for use by the program.")
     
    parser.add_argument("-G", "--genbank_directory", dest="genbank_directory", metavar="DIRECTORY", default='./prokka/',
                help="Folder containing all genbank files for use by the program.")        
                 
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="DIRECTORY", default='./db/',
                help="Folder where the BLAST searchable databases will be stored.")
    
    parser.add_argument("-f", "--filter", dest="filter", metavar="FILE", default='NONE',
                help="File restrictiong which accession numbers this script will process. If no file is provided, filtering is not performed.")
                
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")
                
    parser.add_argument("-p", "--protein", dest="protein", default=True, action='store_false',
                help="Flag to toggle mode from protein to DNA sequences for use as database construction. The default is protein.")
                
    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False,
                help="Suppresses most program text outputs.")
    
    return parser.parse_args()


def check_options(parsed_args):
    '''
    # section of code that checks the infolder entry    
    if os.path.isdir(parsed_args.infolder):
        infolder = parsed_args.infolder
    else:
        print "The folder %s does not exist." % parsed_args.infolder
        sys.exit()
    '''
    # check the genbank folder
    if os.path.isdir(parsed_args.genbank_directory):
        genbank_directory = parsed_args.genbank_directory
    else:
        print "The folder %s does not exist." % parsed_args.genbank_directory
        sys.exit()

    # if the directory that the user specifies does not exist, then the program makes it for them. 
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    outfolder = parsed_args.outfolder
    if outfolder[-1] != '/':
        outfolder = outfolder + '/'
    
    # Check the filter file
    if parsed_args.filter == 'NONE' or os.path.exists(parsed_args.filter):
        filter_file = parsed_args.filter
    else:
        print "The file %s does not exist." % parsed_args.filter
        sys.exit()
    
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)
    
    do_protein = parsed_args.protein
    
    quiet = parsed_args.quiet
      
    #return infolder, outfolder, filter_file, num_proc, do_protein
    return genbank_directory, outfolder, filter_file, num_proc, do_protein, quiet
 
        
#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname) and '.gbk' in fname:
                result.append(fname)
    return result
    
    
#################################################################
#### so yeah.... not gonna do this right now...             #####
#### I feel there is no reason to, but keeping              #####
#### the code in here just in case i am wrong about that    #####
#################################################################

# I need to do an analysis on gc content, and gc skew.
# currently this will return the  organism gc, (mean, SD, varience) of the GC in coding regions
# the results will we retained in a file the calling program opens tab delimited

def GCAnalysis(NC, organism, gc_list, seq, outfile):
    handle = open(outfile, 'a')
    organism_gc = "%3.2f" % GC(seq)
    mean = "%3.2f" % numpy.mean(gc_list)
    #mode = "%5.2f" % numpy.mode(gc_list)
    var = "%5.2f" % numpy.var(gc_list)
    std = "%5.2f" % numpy.std(gc_list)
    handle.write('\t'.join([NC, organism, organism_gc, mean, var, std]) + NEW_LINE)
    

# take the genbank file specified by genbank path, and save the customized result file in the db_directory folder
#def convert_genbank(genbank_path, db_directory, error_fname): #, gc_outfile = 'gc_analysis.txt'):
def convert_genbank(genbank_tuple):
    genbank_path, db_directory, error_fname, do_protein = genbank_tuple
    record_list = []
    seq_record = SeqIO.parse(open(genbank_path), "genbank").next()
    print (seq_record.annotations)
    accession = seq_record.id
    organism = seq_record.annotations['organism'].replace(' ', '_')
    err_log = []
    gc_list = [] # no need for this right now, but leaving in
    # loop over the genbank file
    for fnum, feature in enumerate(seq_record.features):
        err_flag = False
        error_in_field = False
        if feature.type == 'CDS':
            #print dir(feature.location)
            try:
                start = int(feature.location.start)
                stop = int(feature.location.end)
            except:
                error_in_field = True
            
            strand = feature.strand
            dna_seq = seq_record.seq[start:stop]
            #print "dna_seq", type(dna_seq), dna_seq
            gc = GC(dna_seq)
            gc_list.append(gc)
            gc = "%3.2f" % gc
            
            try:
                locus = feature.qualifiers['locus_tag'][0]
            except:
                try:
                    locus = feature.qualifiers['gene'][0]
                except:
                    locus = 'error'
                    print "Error in the organism %s with NC # %s" % (organism, accession)
                    err_flag = True
                    err_log.append([organism, accession])

            if do_protein:
                #seq = seq.translate()
                #print type(seq)
                #print feature.qualifiers.keys()
                #seq = dir(feature)
                try:
                    if 'translation' in feature.qualifiers.keys():
                        prot_seq = Seq(''.join(feature.qualifiers['translation']), IUPAC.protein)
                        #print "prot_seq", type(prot_seq), prot_seq
                        
                        if 'gene' in feature.qualifiers:
                            gene = feature.qualifiers['gene'][0]
                            #record_list.append(SeqRecord(prot_seq, id = '|'.join([accession, organism, locus, gene, str(start), str(stop), str(strand), gc]).replace(' ', ''), description = ''))
                            seq_rec_to_store = SeqRecord(prot_seq, id = '|'.join([accession, organism, locus, gene, str(start), str(stop), str(strand), gc]).replace(' ', ''), description = '')
                        else:
                            #record_list.append(SeqRecord(prot_seq, id = '|'.join([accession, organism, locus, 'unknown', str(start), str(stop), str(strand), gc]).replace(' ', ''),description = ''))
                            seq_rec_to_store = SeqRecord(prot_seq, id = '|'.join([accession, organism, locus, 'unknown', str(start), str(stop), str(strand), gc]).replace(' ', ''),description = '')
                            #print prot_seq
                    else:
                        pass
                        #print "This was not a protein sequence"
                except:
                    print "Error in function convert_genbank(genbank_tuple) from the format_db.py script, unhandled error in the genbank parse."
            else:
                # put something in here that will deal with RNA later, if we plan to go that route.
                pass
            if not error_in_field:
                record_list.append(seq_rec_to_store)
            else:
                print "a record was omitted"
            
            '''        
            #print len(seq)
            if len(seq) < 2:
                #pass
                print "len seq", len(seq)
            
            elif do_protein:
                if 'gene' in feature.qualifiers:
                    gene = feature.qualifiers['gene'][0]
                    record_list.append(SeqRecord(seq, id = '|'.join([accession, organism, locus, gene, str(start), str(stop), str(strand), gc]).replace(' ', ''),
                       description = ''))
                else:
                    record_list.append( 
                      SeqRecord(seq, id = '|'.join([accession, organism, locus, 'unknown', str(start), str(stop), str(strand), gc]).replace(' ', ''),
                       description = ''))
            
            
            else:
                if 'gene' in feature.qualifiers:
                    gene = feature.qualifiers['gene'][0]
                    record_list.append(SeqRecord(seq, id = '|'.join([accession, organism, locus, gene, str(start), str(stop), str(strand), gc]).replace(' ', ''),
                       description = ''))
                else:
                    record_list.append( 
                      SeqRecord(seq, id = '|'.join([accession, organism, locus, 'unknown', str(start), str(stop), str(strand), gc]).replace(' ', ''),
                       description = ''))
                       '''
    #if os.path.isfile(gc_outfile):
    #    os.remove(gc_outfile)
    #GCAnalysis(accession, organism, gc_list, seq_record.seq, gc_outfile)    
    handle = open(error_fname, 'a')
    for i in err_log:
        handle.write('\t'.join(i) + '\n')
        handle.close()
    if not err_flag:
        outpath = db_directory + os.path.splitext(os.path.basename(genbank_path))[0] + '.ffc'
        #print outpath
        out_handle = open(outpath,"w")
        SeqIO.write(record_list, out_handle, "fasta")
        out_handle.close()
        
    if do_protein:
        cmd = "makeblastdb -in %s -dbtype prot" % (outpath)
        #print "got here"
    else:    
        cmd = "makeblastdb -in %s -dbtype prot" % (outpath)
    os.system(cmd)
    #print "Passed main loop"
  
    return outpath, err_flag


def parallel_convert_genbank(file_list, outfolder, num_proc, do_protein, error_fname = "./error_log.txt"):
    
    # Make sure that we have a new error log each time the program is run
    if os.path.isfile(error_fname):
        os.remove(error_fname)
    
    # Package the variables for the convert_genbank function so everything can be run in parallel    
    tuple_list = [(i, outfolder, error_fname, do_protein) for i in file_list]
    
    pool = Pool(processes = num_proc)
    result = dict(pool.map(convert_genbank, tuple_list))


def main():
    
    start = time.time()

    parsed_args = parser_code()
    
    #infolder, outfolder, filter_file, num_proc, do_protein = check_options(parsed_args)
    genbank_directory, outfolder, filter_file, num_proc, do_protein, quiet = check_options(parsed_args)
    
    #flist = returnRecursiveDirFiles(infolder)
    flist = returnRecursiveDirFiles(genbank_directory)

    if filter_file != 'NONE':
        filter_list = [i.strip() for i in open(filter_file).readlines()]
        file_list = [i for i in flist if i.split('/')[-1].split('.')[0] in filter_list]
    else:
        file_list = flist 
    
    #print "do_protein", do_protein
    parallel_convert_genbank(file_list, outfolder, num_proc, do_protein)

    if not quiet: 
        print time.time() - start

    # A successful command could look like this:
    # ./format_db.py -f ./phylo_order.txt
    # ./format_db.py -i /home/dave/Desktop/all_genbank -o ./db1/ -f ./phylo_order.txt
    
if __name__ == '__main__':
    main()
