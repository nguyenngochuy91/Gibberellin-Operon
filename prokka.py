#!/usr/bin/python

import argparse
import os
from shutil import copyfile
def parser_code():

    parser = argparse.ArgumentParser(description='')

    parser.add_argument("-i", "--input_dir", default='./family_genome_folder/',
                help="Folder fasta file (family_genome_folder).")
    parser.add_argument("-o", "--output_dir", default='./prokka/',
                help="Folder that contains gbk file(prooka).")
                
    return parser.parse_args()

#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result
    
if __name__ == '__main__':
    args = parser_code()
    input_dir = args.input_dir
    output_dir = args.output_dir
    result = returnRecursiveDirFiles(input_dir)
    for file_name in result:
        species_name = file_name.split('/')[2]
        print (species_name)
        # create a dir with this name in the prokka dir
        mydir = output_dir+species_name
        cmd ="prokka --outdir "+mydir+" --prefix "+species_name+" "+file_name
        print (cmd)
        os.system(cmd)
        print ("Finish with"+species_name)
#     get the reference gbk file
    reference_dir = output_dir+"reference/"
    print (reference_dir)
    try:
        os.mkdir(reference_dir)
    except:
        print ("reference dir already made")
    copyfile("reference.gbk", reference_dir+"reference.gbk")