#!/usr/bin/env python

from multiprocessing import Pool
import time
import os
import sys
import argparse
import math
from homolog4 import *
from collections import defaultdict
import itertools
from collections import Counter
from Levenshtein import distance
import cPickle as pickle

# Copyright(C) 2015 David Ream, modified by Huy Nguyen
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description="This program will calculate the number of events that an organismal pair do not have in common.  The return is a pickled dictionary.")
    
    parser.add_argument("-i", "--gene_block_file", dest="gene_block_file", default='./regulonDB/gene_block_names_and_genes.txt', metavar="FILE",
                help="The parsed gene block file. The default is the parsed regulonDB operon file.")
                
    parser.add_argument("-I", "--infolder", dest="infolder", default='./optimized_gene_block_1/', metavar="DIRECTORY",
                help="A folder that contains the final gene block results. Gene blocks are grouped by organism, ordered by start, and have spurious BLAST results removed.")
    
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="DIRECTORY", default='./gene_block_distance_matrices/',
                help="Folder where the results will be stored. Default is the folder './gene_block_distance_matrices/'.")

    parser.add_argument("-F", "--gene_block_filter", dest="gene_block_filter", default='NONE', metavar="FILE",
                help="A file that contains the gene blocks that are under investigation.  All others will be omitted from analysis an results.") 
                
    # this is where i will/would add the option to filter organisms.       
    
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")
                
    parser.add_argument("-g", "--max_gap", dest="max_gap", metavar="INT", default = 500, type=int,
                help="Size in nucleotides of the maximum gap allowed between genes to be considered neighboring. The default is 500.")
                            
    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False,
                help="Suppresses most program text outputs.")
    
    return parser.parse_args()


def check_options(parsed_args):

    if os.path.exists(parsed_args.gene_block_file):
        gene_block_file = parsed_args.gene_block_file
    elif parsed_args.gene_block_file == 'NONE':
        gene_block_file = 'NONE'
    else:
        print "The file %s does not exist." % parsed_args.gene_block_file
        sys.exit()

    if os.path.isdir(parsed_args.infolder):
        infolder = parsed_args.infolder
    else:
        print "The folder %s does not exist." % parsed_args.infolder
        sys.exit()
    
    # if the directory that the user specifies does not exist, then the program makes it for them. 
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    outfolder = parsed_args.outfolder
    if outfolder[-1] != '/':
        outfolder = outfolder + '/'
     
    if os.path.exists(parsed_args.gene_block_filter):
        gene_block_filter = parsed_args.gene_block_filter
    elif parsed_args.gene_block_filter == 'NONE':
        gene_block_filter = parsed_args.gene_block_filter
    else:
        print "The file %s does not exist." % parsed_args.gene_block_filter
        sys.exit()
        
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)
    
    # validate the input for the maximum allowed gap
    try:    
        max_gap = int(parsed_args.max_gap)
        if max_gap <= 0:
           print "The gap that you entered %s is a negative number, please enter a positive integer." % parsed_args.max_gap
           sys.exit()
        else:
           pass
    except:
        print "The gap that you entered %s is not an integer, please enter a positive integer." % parsed_args.max_gap
        sys.exit()
        
    quiet = parsed_args.quiet
    
    return gene_block_file, infolder, outfolder, gene_block_filter, num_proc, max_gap, quiet

#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def return_recursive_dir_files(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# This function returns the list of gene block files, filtered by gene block if a filter file is supplied by the user.    
def return_file_list(infolder, filter_file):
    if filter_file == 'NONE':
        return return_recursive_dir_files(infolder)   
    else:
        filter_list = [i.strip() for i in open(filter_file)]
        return [i for i in return_recursive_dir_files(infolder) if os.path.basename(i).split('.')[0] in filter_list]

# This function is no longer used it would seem.    
# This function will take the folder of the gene block result files, and it will return a dictionary that has as a primary key
# the gene block name, the secondary key will be nc. The values of the attributes will then follow.
def parse_operon_result_files(in_folder, dest_folder, gene_block_filter_file):
    file_list = return_file_list(in_folder, gene_block_filter_file)
    
    print "file_list", file_list
    
    result = {}
    #distmat_dict = {}
    
    iddo_result = {}
    
    #gene_dict = make_gene_block_gene_string_dict(gene_block_file)
    #print len(gene_dict)
    
    #print "in_folder", in_folder

    #for i in open(distmat_file).readlines():
    #    nc = i.strip().split('\t')[0]
    #    val = i.strip().split('\t')[1].strip() # fixes a wonky error in distmat.... ugh
    #    distmat_dict.update({nc: val})
    for fname in file_list:
        gene_block = fname.split('/')[-1].split('.')[0]
        #fname = in_folder + f
        print "fname", fname, gene_block
        result.update({gene_block: {}})
        
        iddo_result.update({gene_block: {}})
        
        file_output = []
        summary_info = [i.strip().split('\t')[1:] for i in open(fname).readlines() if i[:2] == '##']
        homolog_entries = []
        
        #print "summary_info", summary_info
        # This section of code has a single purpose. ad-hoc the goodness of rearrangements and distance
        # between grouped genes.
        
        summmary_info = []
        tmp_hlog_list_for_grouping = []
        ignore_list = ['#', '$', '@', '+']
        #for i in [i.strip() for i in open(fname).readlines() if len(i) > 1]: # ugh, corrects for some garbage in the files (occurs once it seems)
        for i in [i.strip() for i in open(fname).readlines() if i.split('\t')[0] == '##']:
        #for i in [i.strip() for i in open(fname).readlines() if len(i) > 1 and i[0] == '#']:
            #if i[:2] == '##':
            if len(i) < 2:
                print "Error in the function parse_operon_result_files from the make_event_distance_matrix script", i, fname
            if i[0] == '#':
                comprehensive_list, group_list = group_homologs(tmp_hlog_list_for_grouping, INTERGENIC_MAX_LENGTH)
                for group in group_list:
                    #print gene_dict[gene_block]['reference_string']
                    rearrangements = return_operon_string_distance(gene_dict[gene_block]['reference_string'], return_group_str(group, gene_block, gene_dict))
                    print gene_dict[gene_block]['reference_string'], return_group_str(group, gene_block, gene_dict), return_gene_block_string_distance(gene_dict[gene_block]['reference_string'], return_group_str(group, gene_block, gene_dict))
                    print "vals recorded", gene_dict[gene_block]['reference_string'], return_group_str(group, gene_block, gene_dict), return_gene_block_string_distance(gene_dict[gene_block]['reference_string'], return_group_str(group, gene_block, gene_dict)), rearrangements
                    print "i line", i, len(i.split('\t'))
                    try:
                        a,nc,c,d,e,f,g,h,i,j,k = i.split('\t') # only interested in a few of the fields here: 
                        distmat_val = distmat_dict[nc]
                        print c
                        common = '_'.join(c.split('_')[:2])
                        print "common", common
                    except:
                        print "Error in line", i, fname
                    
                    
                #print len(tmp_hlog_list_for_grouping), tmp_hlog_list_for_grouping, len(group_list)
                tmp_hlog_list_for_grouping = []
            elif i[0] not in ignore_list:
                tmp_hlog_list_for_grouping.append(return_homolog(i))
                pass
                #print i.split('\t')
        
        
        
        ignore_list = ['#', '$', '@', '+']
        
        #for line in [i.strip() for i in open(fname).readlines() if i[0] not in ignore_list ]:
        #    print line 
        #print homolog_entries
        #print summary_info
        #print "summary_info2", summary_info
        for item in summary_info:
            nc = item[0]
            common = ' '.join(item[1].split('_')[:2])
            if nc in distmat_dict:
                distmat_val = distmat_dict[nc]
                #print nc, common
                line = [nc, common, distmat_val] + [i.strip() for i in item[2:]]
                #print line
                file_output.append(line)
        header = 'NC_Number,Common,Distance,Splitting,Migration,Duplicates,Deletions,Splits,Longest_Group,Total_Fusions,Group_Fusions\n'
        handle = open(dest_folder + gene_block + '.csv', 'w')
        handle.write(header)
        handle.write('\n'.join([','.join(i) for i in file_output]))
        handle.close()
    #print result.keys()    

'''    
#####################################
# No longer used					#
#####################################
# This function will return information about the gene blocks from a parse gene block file.
# The default is the parsed operon file form regulonDB, but the user could also provide one.
# The format is operon_name gene1 gene2 gene3 gene4 etc...  tab delineated.
def return_gene_block_dictionary(flist, gene_block_file):

    # Step 1: Determine information about the gene block(s) that we are investigating
    gene_block_length_dict = {}
    for line in [i.strip().split('\t') for i in open(gene_block_file).readlines()]:
        gene_block = line[0]
        gene_list = line[1:]
        number_genes = len(gene_list)
        gene_block_length_dict.update({gene_block:{'number_genes':number_genes, 'gene_list':gene_list}})
        
    # testing step 1
    print sorted(gene_block_length_dict.keys())
    
    return gene_block_length_dict
'''

# This function will return the information about each organism in a given gene block
# The information will be the accession, and a gene count of each of the gene recovered and the number of groups.
# The returned dict will be keyed {gene_block:{accession:[gene1, gene2, gene3, etc...]}.  This will be used
# in later steps to generate the gene and split counts used to create the distance matrix in later steps.
def return_gene_block_organism_data(flist):
    result = {}
    
    for fname in flist:
        # this will store the homologs in a dict {accession:[hlog1, hlog2, etc]} temporarily. 
        org_dict = {}
        gene_block = os.path.basename(fname).split('.')[0]
        #print "gene_block", gene_block
        
        for line in [i.strip() for i in open(fname).readlines()]:
            # catch some errors if they occur and report a useful message as a result.
            try:
                hlog = Homolog.from_blast(line)
            except:
                print "Error in function return_gene_block_organism_data from script make_event_distance_matrix.py", line
            try:
                accession = hlog.accession()
            except:
                print "There was an error in the homolog class in function return_gene_block_organism_data from script make_event_distance_matrix.py"
            predicted_gene = hlog.blast_annotation()
            
            # store the homolog where it belongs
            if accession in org_dict.keys():
                org_dict[accession].append(hlog)
            else:
                org_dict.update({accession:[hlog]})

        # Now make sure that everything is put neatly in order by start position.
        for accession in sorted(org_dict.keys()):
            h_list = org_dict.pop(accession)
            h_list.sort(key=lambda x: x.start())
            org_dict.update({accession:h_list})
        
        # Store the dict as an entry in the result.  
        result.update({gene_block:org_dict})

    return result

# adding a list of keys that are not the names of genes. These keys exist to keep track of non-gene information about the 
# gene blocks that were are investigating.
IGNORE_LIST = ['groups']


# This function strips the necessary information about the homologs that we report per gene block/organism pair
# so that event determination between different organisms can be performed.  No events are determined in this stage,
# simply information gathering, for later compairision.
# If we wish to add new events to our method, this is the first place that we will need to add new code.
def make_event_reporting_data_structure(gene_block_dict, max_gap):
    result = {}
    
    for gene_block in sorted(gene_block_dict.keys()):
        #print gene_block
        result.update({gene_block:{}})
        for organism in gene_block_dict[gene_block].keys():
            #print organism
            # we use this list twice, so storing it
            list_homologs = gene_block_dict[gene_block][organism]
            # returns a count of each gene that was in the homolog list as a dictionary
            org_gene_cnt = dict(Counter([i.blast_annotation() for i in  list_homologs]))
            # returns a list of gene blocks, the number of which can be used to determine splits later. 
            # the special name 'groups' is how that can be determined.
            gene_block_groupings = homolog_list_grouping_function(list_homologs, max_gap)
            groups = len(gene_block_groupings)
            # save the number of groups to the result dict for the organism under investigation
            org_gene_cnt.update({'groups':groups})
            # Store the result
            result[gene_block].update({organism:org_gene_cnt})
    return result
                

# This function will return the number of splits that two organisms do not share.
def return_splits(data_struct, org1, org2):
    result = 0
    if org1 != org2:
        org1_groups = data_struct[org1]['groups']
        org2_groups = data_struct[org2]['groups']
        result = int(math.fabs(org1_groups - org2_groups))
    else:
        pass
        
    return result

# this function will return the number of deletions as the number of unique genes they do not share
def return_deletions(data_struct, org1, org2):
    if org1 != org2:
        # Removing the key 'genes', and any other non gene information in the reporting data_struct using the ignore_list
        #gene_list_org1 = data_struct[org1].keys()
        gene_list_org1 = [i for i in data_struct[org1].keys() if i not in IGNORE_LIST]
        #gene_list_org2 = data_struct[org2].keys()
        gene_list_org2 = [i for i in data_struct[org2].keys() if i not in IGNORE_LIST]
        #unique_gene_list = [i for i in list(set(gene_list_org1) - set(gene_list_org2)) if i not in IGNORE_LIST]
        unique_gene_list = list(set(gene_list_org1 + gene_list_org2))
        
        # Here we determine how many deletions there are in each organism. The total of these deletions is the number
        # of deletions between the organisms.
        deletions_in_org1 = len(list(set(unique_gene_list)- set(gene_list_org1)))
        deletions_in_org2 = len(list(set(unique_gene_list)- set(gene_list_org2)))
        return deletions_in_org1 + deletions_in_org2
        
    else:
        return 0
        
# this function will return the number of deletions as the number of duplicated genes
def return_duplications(data_struct, org1, org2):
    result = 0
    if org1 != org2:
        gene_list_org1 = data_struct[org1].keys()
        gene_list_org2 = data_struct[org2].keys()
        unique_gene_list = [i for i in list(set(gene_list_org1) - set(gene_list_org2)) if i not in IGNORE_LIST]
        for gene in unique_gene_list:
            duplications = 0
            if gene in data_struct[org1].keys():
                gene1_copy_number = data_struct[org1][gene]
            else:
                gene1_copy_number = 0
            if gene in data_struct[org2].keys():
                gene2_copy_number = data_struct[org2][gene]
            else:
                gene2_copy_number = 0
                
            #########################################################################################################################
            # This is where we calculate duplications.  The issues here are that we first consider a deletion BEFORE a duplication.	# 
            #  Therefore if organism 1 lacks a gene, for this explanation call it gene 'A', and organism 2 has two copies of it, we	#
            # will say that the first copy of gene A in organism 1 is deleted with respect to organism 2.  Similarily the second 	#
            # of this gene is a duplication.  This will yield the situation where we count one deletion and one duplication for this#
            # situation.																											#
            #########################################################################################################################
            
            gene_list = [gene1_copy_number, gene2_copy_number]
            
            # if one organisms is missing the gene entirely
            if min(gene_list) == 0:
                duplications =  sum(gene_list) - 1
            # if both organisms have the gene
            else:
                duplications =  int(math.fabs(gene1_copy_number - gene2_copy_number))
            
            result += duplications   
                    
    else:
        # org 1 and org2 are the same, so we do nothing
        pass
        
    return result
 

# This function will return an all vs. all pickled dict of the format: 
# {gene_block:{NC1:{NC2:{event1:numeric, event2:numeric, etc:numeric}}}}
def return_event_counts_as_dict(event_reporting_data_structure):
    result = {}
    
    for gene_block in sorted(event_reporting_data_structure.keys()):
        #print "gene_block ", gene_block 
        result.update({gene_block:{}})
        org_list = event_reporting_data_structure[gene_block]
        # use the magic that is itertools to make an iterable for all combinations of two organisms for compairison 
        for pair in itertools.combinations_with_replacement(org_list, 2):
            #print "Pair", pair
            org1, org2 = pair
            deletions = return_deletions(event_reporting_data_structure[gene_block], org1, org2)
            duplications = return_duplications(event_reporting_data_structure[gene_block], org1, org2)
            splits = return_splits(event_reporting_data_structure[gene_block], org1, org2)
            
            #print "org1", org1, "org2", org2 
            
            #print "Deletions", deletions, "Duplications", duplications, "Splits", splits
            
            # There is a faster implementation that uses a try except... which i have been kind enough to include if speed is you thing.
            # on small lists of organisms it really does not matter
            '''if org1 in result[gene_block].keys():
                result[gene_block][org1].update({org2:{'deletions':deletions, 'duplications':duplications, 'splits':splits}})
            else:
                result[gene_block].update({org1:{org2:{'deletions':deletions, 'duplications':duplications, 'splits':splits}}})
            '''
            if org1 in result[gene_block].keys():
                result[gene_block][org1].update({org2:{'deletions':deletions, 'duplications':duplications, 'splits':splits}})
            else:
                result[gene_block].update({org1:{org2:{'deletions':deletions, 'duplications':duplications, 'splits':splits}}})
            
            
            '''  
            try:
                result[gene_block][org1].update({org2:{'deletions':deletions, 'duplications':duplications, 'splits':splits}})
            except:
                result[gene_block].update({org1:{org2:{'deletions':deletions, 'duplications':duplications, 'splits':splits}}})
            '''
            
    return result

# This function will return the rearrangement distance between two organisms by using the levinstein edit distance
# as a rough metric.  
def return_rearrangements(sorted_gene_block_results, org1, org2):
    if org1 != org2:
        pass
    else:
        return 0 
               
# Take a list of unique genes between two organisms and return a mapping for gene name to a single character.  
# This will be used to determine the Levenshtein edit distance. 
def make_gene_string_dict(gene_list):
    result = {}
    for gene, index in zip(gene_list, range(0,len(gene_list))):
        #print gene
        result.update({gene:chr(65+index)})
    #print operon, result[operon], len(result)
    return result

# This function will take a list of ordered homologs, and groups them by a max_gap constraint.
# The return is a list of lists. Single genes and gene blocks will both be represented as groups.
def homolog_list_grouping_function(list_homologs, max_gap):
    result = []
    neighborhood = [list_homologs[0]]
    
    for i in list_homologs[1:]:
        #look at current
        start = neighborhood[-1].start() #start = list_homologs[i].start()
        stop = neighborhood[-1].stop() #stop = list_homologs[i].stop()
        # look at next
        start_n = i.start() #start_n = list_homologs[i+1].start()
        stop_n = i.stop() #stop_n = list_homologs[i+1].stop()
        
        # We have found neighboring genes, as defined by max_gap
        if math.fabs(start - stop_n) < max_gap or math.fabs(stop - start_n) < max_gap:
            neighborhood.append(i)
        # These genes do not neighbor eachother
        else: 
            result.append(neighborhood)
            neighborhood = [i]
    result.append(neighborhood)
    #print list_homologs[0].organism(), "result", result, "neighborhood_found ", neighborhood_found  
    return result
    
def main():
    start = time.time()
    
    parsed_args = parser_code()
    
    gene_block_file, infolder, outfolder, gene_block_filter, num_proc, max_gap, quiet = check_options(parsed_args)
    
    if not quiet:
        print gene_block_file, infolder, outfolder, gene_block_filter, num_proc, max_gap
    
    # Step 1: Filter the list of gene blocks from the input folder that we are interested in, 
    # and return information about them for use in later steps.
    file_list = return_file_list(infolder, gene_block_filter)
    
    # Step 2: Parse the results using the homolog class, and store result as a dictionary that is keyed on the accession number.
    # we use the accession number since it is the organisimal unique that the rest of the scripts in the project use.
    
    sorted_gene_block_results = return_gene_block_organism_data(file_list)
    
    # Step 3: Calculate the gene counts and number of splits. 
    # In later versions of the script, it may be necessary to add more events, and this code will have to adapt to that.
    # This is the point of insertion for new code to carry the raw information that is used later for event determination.
    
    # TODO: make commandline argument!
    #max_gap = 500
    
    # return a dict that contains the raw information needed to make a between organismal compairison of the events used in our method.
    event_reporting_data_structure = make_event_reporting_data_structure(sorted_gene_block_results, max_gap)
    

    # Step 4: Actually do the compairison and report the result.
    # It will be keyed {gene_block:{NC1:{NC2:{event1:numeric, event2:numeric, etc:numeric}}}}
    
    event_count_dict = return_event_counts_as_dict(event_reporting_data_structure)
    
    #print "Length", len(event_count_dict.keys()), event_count_dict.keys()
    
    #outfile = './event_dict.p'
    outfile = outfolder + 'event_dict.p'
    pickle.dump(event_count_dict, open(outfile, 'w'))
    
        
    #parallel_list_param = [(i, outfolder, max_gap, e_val) for i in file_list]
    
    #parse_operon_result_files(infolder, outfolder, operon_filter_file)

    if not quiet:
        print time.time() - start


if __name__ == '__main__':
    main()   
