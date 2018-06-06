#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : blasting each 
    Start   : 09/28/2016
    End     : /2016
'''
import os
import argparse
def parser_code():

    parser = argparse.ArgumentParser()                  
    parser.add_argument("-o", "--gene_blast", dest="gene_blast", metavar="DIRECTORY",
                help="the blast info of gene blast for cyp114 and cyp115, Fd,SDR,GGPS in fa format")    
                            
    return parser.parse_args()

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

# function to parse and only keep those with percent indentity greater or equal to 70%
def filter_blast(myfile):
    infile = open(myfile,'r')
    lines = infile.readlines()
    infile.close()
    outfile = open(myfile,'w')
    for line in lines:
        item = line.split('\t')
        locus = item[0].split('|')[2]
        if locus != 'XOC_0083' and float(item[2]) >= 70 or locus == 'XOC_0083':
            outfile.write(line)
    outfile.close()
    
def parse_written(target,query,blast_dir):
    cmd1 = ('blastp -query '+query+' -outfmt 6 -out '+blast_dir+target.split('.ffc')[0]+
        ' -subject db/'+target+' -evalue 1e-10 -num_threads 8')
    os.system(cmd1)
    #print ("cmd1",cmd1)
    print ('Finish blasting:',target)

#    filter_blast(blast_dir+target)    
#    print ('Finish filter blast:',target)
    
    

def main():
    args  = parser_code()
    try:
        os.mkdir("gene_blast")
    except:
        print "gene_blast dir already created"
    blast_dir  = "gene_blast/"+args.gene_blast+"/"
    try:
        os.mkdir(blast_dir)
    except:
        print blast_dir +" dir already created"
    gene  = args.gene_blast+".fa"
    result = [i for i in return_recursive_dir_files('db') if i.split('/')[-1].split('.')[-1] == 'ffc']
    print len(result)
    for file in result:
        name = file.split('/')[-1]
        parse_written(name,gene,blast_dir)
    #blast_parse()

main()