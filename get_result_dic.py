#!/usr/bin/python

from homolog4 import *

import os
import sys
import argparse
import time
def parser_code():

    parser = argparse.ArgumentParser()
                     
    
    parser.add_argument("--input","-i", default="./optimized_gene_block_1/gibberellin.txt",help="optimized gene block ")
    parser.add_argument("-gb", "--gene_blast", dest="gene_blast", metavar="DIRECTORY", default='./gene_blast/',
                help="the blast info of gene blast for cyp114 and cyp115")    
    parser.add_argument("-a", "--accession", default='accession_to_common.txt',
                help="accesion number to map from accession number to name")  
    parser.add_argument("--output","-o", default="./result_dic/gibberellin",
                help="where the result be stored")
    parser.add_argument("--threshold","-t", default=350,
                help="length threshold to determine whether a pseudo or not, default at 300")
    parser.add_argument("--gene_name","-g",default="gene_block_names_and_genes.txt",
                        help = "gene name of the operon")
    parser.add_argument("-f", "--filter_file", dest="filter_file",  default='potential_fusion.txt',
                help="Filter file, default as potential file, if NONE then not doing the parse filter")
                            
    return parser.parse_args()


# parse accession_to_common.txt
def parse_accession(myfile):
    infile = open(myfile,'r')
    result = {}
    for line in infile.readlines():
        line = line.strip().split(',')
        if line[0]!= 'CP003057.2':
            result[line[0]] = line[1]
        else:
            result[line[0]] = 'reference'
    return result
#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result
    
# parsing the filter file, here is the potential_fusion.txt, get the species name out
# input: filter file
# output: list of species that have potential fusion
def parse_filter(filter_file):
    result = []
    infile = open(filter_file,'r')
    for line in infile.readlines():
        line = line.strip()
        info = line.split(':')[1]
        info = info.split(' ')[0]
        result.append(info)
    return result

def parse_GGPS2(gene_blast_ggps):
    result =  [i.split('/')[-1] for i in  returnRecursiveDirFiles(gene_blast_ggps)]

    dic ={}
    for species in result:
        handle = open(gene_blast_ggps+species,"r")
        for line in handle.readlines():
            line             = line.split('\t')
            target           = line[1]
            percent_indentity = line[2]
            info             = target.split("|")
            startPos         = info[4]
            endPos           = info[5]
            strandDirection  = info[6]
            # get the specie name
            target = info[0]
            if target in dic:
                if float(dic[target][0])<float(percent_indentity):
                    dic[target] = [percent_indentity,startPos,endPos,strandDirection,False]
            else:
                dic[target] = [percent_indentity,startPos,endPos,strandDirection,False]

    return dic
###############################################################################
## dueling with pseudo gene , and categorize cyp114 vs cyp115
# using the result from the filter file, and the 2 dir of cyp114 and cyp115 to get the dictionary
# that store info for each species, which are the poteintial genes
# input: result list, 2 dir paths
# output: dic
def retrieve_info(cyp114_dir, cyp115_dir):
    result =  [i.split('/')[-1] for i in  returnRecursiveDirFiles(cyp114_dir)]

    dic={}
    for species in result:
        dic[species] = {114:{},115:{}}
        infile_114 = open(cyp114_dir+species,'r')
        get_info_into_dic_1(infile_114,dic[species][114])

        infile_115 = open(cyp115_dir+species,'r')
        get_info_into_dic_1(infile_115,dic[species][115])
    # print dic
    return dic

# function that given a handle, read out info into a dictionary
# input  : handle, dic
# output : None 

def get_info_into_dic_1(handle,dic):
    for line in handle.readlines():
        line             = line.split('\t')
        target           = line[1]
        percent_indentiy = line[2]
        alignment_length = line[3]
        mismatch_number  = line[4]
        gap_opening_num  = line[5]
        
        # get the target name
        target   = target.split('|')[2]
        dic[target] = [percent_indentiy,alignment_length,mismatch_number,gap_opening_num]
        
# from the dic, we will find the true pseudo gene for each species
# input  : dic from above
# output : dic that has key is species name, value is the pseudo gene and 1 real gene
    
def get_real_gene(dic,real_gene,fake_gene):
    result = {}
    for species in dic:
        result[species] = []
        for candidates in dic[species][real_gene]:
            info = dic[species][real_gene][candidates]
            if candidates in dic[species][fake_gene]: # check if it is a homolog to cyp114 as well
                # do the comparison
                info_cyp114 = dic[species][fake_gene][candidates]
                info_cyp115 = dic[species][real_gene][candidates]
                if float(info_cyp115[0])- float(info_cyp114[0])> 0:
                    true_length_cyp115 = int(info_cyp115[1])- int(info_cyp115[2]) - int(info_cyp115[3])
                    true_length_cyp114 = int(info_cyp114[1])- int(info_cyp114[2]) - int(info_cyp114[3])
                    if true_length_cyp115 > true_length_cyp114:
                        result[species].append(candidates)
            else:
                result[species].append(candidates)
#    for species in result:
#        print species, result[species]
    return result
    
    
###############################################################################
## dueling with fusion
# using the result from the filter file, and the 3 dir of cyp114 and fd and sdr to get the dictionary
# that store info for each species, which are the poteintial genes
# input: result list, 3 dir paths
# output: dic
def retrieve_info_filter(result, cyp114_dir, fd_dir,sdr_dir):
    dic={}
    for species in result:
        dic[species] = {"114":{},"fd":{},"sdr":{}}

        infile_114 = open(cyp114_dir+species,'r')
        get_info_into_dic(infile_114,dic[species]["114"])
        
        infile_fd = open(fd_dir+species,'r')      
        get_info_into_dic(infile_fd,dic[species]["fd"])
        
        infile_sdr = open(sdr_dir+species,'r')
        get_info_into_dic(infile_sdr,dic[species]["sdr"]) 

    return dic
    
# function that given a handle, read out info into a dictionary
# input  : handle, dic
# output : None 

def get_info_into_dic(handle,dic):
    for line in handle.readlines():
        line             = line.split('\t')
        target           = line[1]
        percent_indentiy = line[2]
        alignment_length = line[3]
        target_start     = line[8]
        target_stop      = line[9]
        
        # get the target name
        target   = target.split('|')[2]
        dic[target] = [percent_indentiy,alignment_length,target_start,target_stop]
        
# given info from fd_info list, and either sdr or cyp114 info, determine wheter
# the blast hits are overlapped, if not then just return as it is
# if not, modify the start and stop point
# input: (fd_info, query_info)
# output: (fd_info, query_info)
def check_overlap(fd_info,query_info):
    fd_start   = int(fd_info[2])
    fd_info[2] = fd_start
    fd_stop    = int(fd_info[3])
    fd_info[3] = fd_stop
    query_start   = int(query_info[2])
    query_info[2] = query_start
    query_stop    = int(query_info[3])
    query_info[3] = query_stop
    if query_start < fd_start: # means that fd is fused into 3'end of the query
        # check if overlap
        if query_stop >= fd_start: # have to change
            middle = (query_stop+fd_start)//2
            query_info[3] = middle
            fd_info[2]    = middle - 1
    else: # means that fd is fused into the 5' end of the 5' end of the query
        if fd_stop >= query_start:
            middle = (fd_stop+query_start)//2
            query_info[2] = middle
            fd_info[3]    = middle - 1
        
    return fd_info,query_info

# from the dic, we will find the fusion if it happens
# input  : dic from above
# output : dic that has key is species name, value is the fusion
    
def modify_prelim_dic(dic):
    for species in list(dic.keys()):
        fd_list = list(dic[species]['fd'].keys())
        if len(fd_list)!=0:
            dic_114 = {}
            sdr_dic = {}
            for locus_name in fd_list: # we assume that we guarantee to find 1 from either sdr or cyp114
                found_sdr = False
                found_114 = False               
                fd_locus_info  = dic[species]["fd"][locus_name]
                if locus_name in dic[species]['sdr']:
                    sdr_info  = dic[species]['sdr'][locus_name]
                    found_sdr = True
                    # check if overlap? if not life is happy
                    fd_locus_info,sdr_info = check_overlap(fd_locus_info,sdr_info)
                    # update the dic_114
                    sdr_dic[locus_name]   = sdr_info

                if locus_name in dic[species]['114']:
                    info_114  = dic[species]['114'][locus_name]
                    found_114 = True
                    # check if overlap? if not life is happy
                    fd_locus_info,info_114 = check_overlap(fd_locus_info,info_114)
                    # update the sdr_dic
                    dic_114[locus_name] = info_114
            
            if found_114:
                dic[species]['114'] = dic_114
            else:
                del dic[species]['114']
            if found_sdr:
                
                dic[species]['sdr'] = sdr_dic
            else:
                del dic[species]['sdr']
        else:
            del dic[species]
            
    return dic

def parse(gene_name):
    result = {}
    alphabet= ["k","a","b","c","d","e","f","g","h","i","j"]
    infile = open(gene_name,'r')
    line = infile.readline().strip()
    line = line.split('\t')[1:]
    for i in range(len(line)):
        result[line[i]] = alphabet[i]
    return result
    
    

def format_operon(operon,threshold):
    result = {}
    infile = open(operon,'r')
    for line in [i.strip() for i in infile.readlines()]:
        hlog = Homolog.from_blast(line)
        accession = hlog.accession()
        start = hlog.start()
        end = hlog.stop()
        strand = hlog.strand()
        # gene_name = hlog.blast_annotation()
        locus               = hlog.query_locus()
        identity            = hlog.percent_ident()
        aligned_length      = hlog.aligned_length()
        locus_name          = hlog.locus()
        boolean = aligned_length>=threshold
        
        # print "identity",identity
        if accession in result.keys():
            result[accession].append((locus, start, end, strand,identity,boolean,locus_name))
        else:
            result.update({accession:[(locus, start, end, strand,identity,boolean,locus_name)]})

    return result

# function to quickly return string
def to_string(name,startPos,endPos,strandDirection,identity,boolean):
    return name + ','+ str(startPos)+ ','+ str(endPos) + ','+ str(strandDirection)+','+str(identity)+','+str(boolean)+'\t'
    
if __name__ == "__main__":
    args                  = parser_code()
    accesion_file         = args.accession
    filter_file           = args.filter_file
    gene_blast            = args.gene_blast
    gene_blast_cyp114     = gene_blast+'CYP114/'
    gene_blast_cyp115     = gene_blast+'CYP115/'
    gene_blast_fd         = gene_blast+'FD/'
    gene_blast_sdr        = gene_blast+'SDR/'
    gene_blast_ggps       = gene_blast+'GGPS2/'
    GGPS2_dic             = parse_GGPS2(gene_blast_ggps)
    try:
        os.mkdir("result_dic")
    except:
        print "result_dic already created"
    operon = args.input
    output = args.output
    gene_name = args.gene_name
    threshold = args.threshold
    gene_dic = parse(gene_name)
    accesion_dic = parse_accession(accesion_file)
    # print gene_dic
    operon_dic = format_operon(operon,threshold)
#    print operon_dic
    dic_species_to_cyp114_cyp115 = retrieve_info(gene_blast_cyp114,gene_blast_cyp115)
    result_115                   = get_real_gene(dic_species_to_cyp114_cyp115,115,114)
    result_114                   = get_real_gene(dic_species_to_cyp114_cyp115,114,115)
#    print "result_115: ",result_115
#    print "result_114: ",result_114
    potential_fusion_species  = parse_filter(filter_file)
    fd_into_114 = []
    fd_into_sdr = []
    dic_species_to_fusion     = retrieve_info_filter(potential_fusion_species,
                                                     gene_blast_cyp114,
                                                     gene_blast_fd,
                                                     gene_blast_sdr)
    # now filter the dic_species_to_pseudo based on the locus of fd, 
    # remove species from dic if no fd found
    dic_species_to_fusion     = modify_prelim_dic(dic_species_to_fusion)   
    
    #print (dic_species_to_fusion)
    # write out to file result
    fileout = open(output,'w')
    string =""
    for key in sorted(gene_dic):
        string += key+","+gene_dic[key]+'\t'
    string+='\n'
    for species in operon_dic:
        string +=species+':'
        for l in operon_dic[species]:

            geneName=l[0]
            startPos=int(l[1])
            endPos=int(l[2])
            strandDirection=l[3]
            identity = l[4]
            boolean = l[5]
            locus   = l[6]
            name    = accesion_dic[species]
            # using this locus to check for cyp114 and cyp115
            # correct the cyp114 and cyp115 and the fusion
            if name == "reference":
                string += to_string(l[0],startPos,endPos,strandDirection,identity,boolean)
            else:
                
                if l[0] == 'XOC_0083':
                    if locus in result_115[name]: # check if it actually cyp115 instead of cyp114
                        string += to_string('XOC_0085',startPos,endPos,strandDirection,identity,boolean)
                    else: # if it is actually cyp114, check if there is a fusion of Fd into cyp114
                        # if species has no fusion or fusion does not involve cyp114
                        if name not in dic_species_to_fusion or (name in dic_species_to_fusion and "114" not in dic_species_to_fusion[name]):
                            string += to_string(l[0],startPos,endPos,strandDirection,identity,boolean)
                        else:
        
                            info_114 = dic_species_to_fusion[name]["114"]
                            info_fd  = dic_species_to_fusion[name]["fd"]
                            for item in info_fd:
                                # write cyp_114 new info
                                if item == locus:
                                    string += to_string(l[0],startPos+3*(info_114[locus][2]-1),startPos+3*(info_114[locus][3])-1,strandDirection,identity,boolean)
                                    # write fd new info                        
                                    string += to_string('XOC_0082',startPos+3*(info_fd[locus][2]-1),startPos+3*(info_fd[locus][3])-1,strandDirection,identity,boolean)
                elif l[0] == 'XOC_0085':  # checking fusion for fd into cyp114
                    if locus in result_114[name]:# 
                        # if species has no fusion or fusion does not involve cyp114
                        if name not in dic_species_to_fusion or (name in dic_species_to_fusion and "114" not in dic_species_to_fusion[name]):
                            string += to_string('XOC_0083',startPos,endPos,strandDirection,identity,boolean)
                        else:
                            info_114 = dic_species_to_fusion[name]["114"]
                            info_fd  = dic_species_to_fusion[name]["fd"]
                            fd_into_114.append(name)
                            for item in info_fd:
                                if item == locus:
                                    # write cyp_114 new info                       
                                    string += to_string('XOC_0083',startPos+3*(info_114[locus][2]-1),startPos+3*(info_114[locus][3])-1,strandDirection,identity,boolean)
                                    # write fd new info                        
                                    string += to_string('XOC_0082',startPos+3*(info_fd[locus][2]-1),startPos+3*(info_fd[locus][3])-1,strandDirection,identity,boolean)
                    else:
                        string += to_string(l[0],startPos,endPos,strandDirection,identity,boolean)
                elif l[0] == 'XOC_0081': # checking fusion for fd into Sdr
                    if name not in dic_species_to_fusion or (name in dic_species_to_fusion and "sdr" not in dic_species_to_fusion[name]):
                        string += to_string('XOC_0081',startPos,endPos,strandDirection,identity,boolean)
                    else:
                        info_sdr = dic_species_to_fusion[name]["sdr"]
                        info_fd  = dic_species_to_fusion[name]["fd"]
                        fd_into_sdr.append(name)
                        for item in info_fd:
                            if item == locus:
                                # write cyp_114 new info
                                string += to_string('XOC_0081',startPos+3*(info_sdr[locus][2]-1),startPos+3*(info_sdr[locus][3])-1,strandDirection,identity,boolean)
                                # write fd new info                        
                                string += to_string('XOC_0082',startPos+3*(info_fd[locus][2]-1),startPos+3*(info_fd[locus][3])-1,strandDirection,identity,boolean)
                else:
                    string += to_string(l[0],startPos,endPos,strandDirection,identity,boolean)
        
        # check if this specie has the XOC_0086 which is GGPS2
        if species in GGPS2_dic:
            percent_identity = float(GGPS2_dic[species][0])
            startPos = GGPS2_dic[species][1]
            endPos   =GGPS2_dic[species][2]
            strandDirection = GGPS2_dic[species][3]
            if percent_identity >= 70:
                string +=to_string('XOC_0086',startPos,endPos,strandDirection,percent_identity,True)
        string +='\n'
    fileout.write(string)
    fileout.close()
    outfile = open("fusion_categorization.txt","w")
    outfile.write("Species might have gene Fd fused to gene CYP114:\n")
    for item in fd_into_114:
        outfile.write(item+"\n")
    outfile.write("\n")
    outfile.write("Species might have gene Fd fused to gene Sdr:\n")
    for item in fd_into_sdr:
        outfile.write(item+"\n")   
    outfile.close()