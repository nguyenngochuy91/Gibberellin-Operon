#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : For csv file, get the genbank data file
    Start   : 09/28/2016
    End     : /2016
'''

from Bio import Entrez
from bs4 import BeautifulSoup # for parsing website 
import urllib 
import os
import time 
directory = './family_contigs_folder/'
try:
    os.mkdir(directory)
except:
    print ("family_contigs_folder is already created")
web = 'https://www.ncbi.nlm.nih.gov/genome/?term='
assembly = 'https://www.ncbi.nlm.nih.gov/gquery/?term='
def read_name(myfile):
    infile = open(myfile,'r')
    lines = infile.readlines()
    names =[]
    for line in lines:
        names.append(line[:-1])
    return names
start = time.time()
names = read_name('name')
    
# getting ftp from assembly database
def query_assembly(myfile,name):
    name_dic = {}
    for name in names:
        sub_info = name.split(' ')
        organism = ' '.join(sub_info[:2])
        strain   = ' '.join(sub_info[2:])
        name_dic[name] =[organism,strain]
    my_dic ={}
    infile = open(myfile,'r')
    for line in infile.readlines()[2:]:
        info = line.split('\t')
        for name in names:
            if name == info[7]:
                my_dic[name] = info[19]+'/'+info[19].split('/')[-1]+'_genomic.fna.gz'
                names.remove(name)
            elif name_dic[name][0] in info[7] and name_dic[name][1] in info[8]:
                my_dic[name] = info[19]+'/'+info[19].split('/')[-1]+'_genomic.fna.gz'
                names.remove(name)
            elif info[7] in name and name_dic[name][1].split(' ')[-1] in info[8]:
                my_dic[name] = info[19]+'/'+info[19].split('/')[-1]+'_genomic.fna.gz'
                names.remove(name)
    return my_dic        
# running on assembly sumary
my_dic = query_assembly('assembly_summary.txt',names)
bioproject =[]
for name in names:
    if name not in my_dic:
        bioproject.append(name)  
        

# download our file
def download(my_dic):
    for name in my_dic:
        # download command into the file with right name
        file_name = directory+'_'.join(name.split(' '))+'.gz'
        cmd1 = 'curl -o '+ file_name+ ' '+ my_dic[name]
        print (cmd1)
        os.system(cmd1)
        # unzip the file
        cmd2 = 'gunzip '+file_name
        print (cmd2)
        os.system(cmd2)
        print ("Successfully download assembly of ",name)
# actually download the data 
download(my_dic)
stop = time.time()
print ("Quite long time eh: ",stop - start)