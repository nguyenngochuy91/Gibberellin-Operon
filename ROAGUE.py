#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : the whole pipeline to generate the result
    Start   : 09/28/2017
    End     : /2017
'''
import argparse
import os
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--contigs","-c",help = "contigs length (default as 10000)",default = 10000)
    parser.add_argument("--model","-m",help = "word model file")
    parser.add_argument("--wlcut","-l", help = "the maximum number of words in any superword (a positive integer)")
    parser.add_argument("--pseudo","-p", help = "include pseudo of gene CYP115 in reconstruction or not(y,n)",default = "y")
    args = parser.parse_args()
    return args
    
def main():
    args      = get_arguments()
    includePseudo = args.pseudo
    contigs = args.contigs
    cmd1 = "./get_fasta.py"
    cmd2 = "./combine_contigs.py -t "+str(contigs)
    cmd3 = "./prokka.py"
    cmd4 = "./format_db.py"
    cmd5 = "./blast_script_1.py"
    cmd6 = "./blast_parse.py"
    cmd7 = "./filter_operon_blast_results.py"
    cmd8 = ["./blast_1_gene.py -o CYP114","./blast_1_gene.py -o CYP115","./blast_1_gene.py -o FD",
            "./blast_1_gene.py -o SDR","./blast_1_gene.py -o GGPS2"]
    cmd9 = "./get_result_dic.py "
    cmd10 = "./ancestral_reconstruction/convert.py -i result_dic/ -o ancestral_reconstruction/new_result/ -l XOC_0085"
    cmd11 = "./create_newick_tree.py"
    cmd12 = "./create_operon_tree.py "
#    # cmd1
#    print ("Executing command: {} \n {} \n".format(cmd1," find all the assembly file of such species using the assembly_summary and download the assembly file from ncbi"))
#    os.system(cmd1)
#    
##     cmd2
#    print ("Executing command: {} \n {} \n".format(cmd2,"  combine contigs given threshold "+str(contigs)))
#    os.system(cmd2)
#    
##     cmd3
#    try:
#        os.mkdir("prokka")
#    except:
#        print ("dir prokka already exists")
#    print ("Executing command: {} \n {} \n".format(cmd3,"   using prokka to annotate "))
#    os.system(cmd3)
#  
#     cmd4
#    try:
#        os.mkdir("db")
#    except:
#        print ("dir db already exists")
#    print ("Executing command: {} \n {}\n".format(cmd4,"  Build a db of the cds for each species to blast it "))
#    os.system(cmd4)
#    
##     cmd5
#    try:
#        os.mkdir("blast_result_1")
#    except:
#        print ("dir blast_result_1 already exists")
#    print ("Executing command: {} \n {}\n".format(cmd5,"  blast file in ./db/ vs gene_block_query.fa"))
#    os.system(cmd5)
#        
#    # cmd6
#    try:
#        os.mkdir("blast_parse_1")
#    except:
#        print ("dir blast_parse_1 already exists")
#    print ("Executing command: {} \n {}\n".format(cmd6,"  blast parse"))
#    os.system(cmd6)
#        
#        
##      cmd7 
#    print ("Executing command: {} \n {}\n".format(cmd7,"   filter to get the best operon"))
#    os.system(cmd7)
#    # cmd8
#    print ("create gene blast for CYP114, CYP115, FD, SDR, and GGPS2 \n")
#    for item in cmd8:
#        print ("Executing command: {} \n".format(item))
#        os.system(item)
#    
##     cmd9
#    print ("Executing command: {} \n {}\n".format(cmd9,"  create the result_dic file for the operon"))
#    os.system(cmd9)
#    
#    #cmd10
#    print ("Executing command: {} \n {}\n".format(cmd10," this take in the gibberellin file and check for the pseudo gene and set it to p (naively using only the percentage)"))
#    os.system(cmd10)
#    
#    #cmd11
#    print ("Executing command: {} \n {}\n".format(cmd11," generate species tree (using gene rpoB)"))
#    os.system(cmd11)
#    #cmd12
#    print ("Executing command: {} \n {}\n".format(cmd12," generate operon tree (concatenate all gnees)"))
#    os.system(cmd12)
#    
    
    ######
    ## ancestral reconstruction

    if includePseudo in "Yy":
        cmd13 = "./ancestral_reconstruction/reconstruction.py -i ancestral_reconstruction/new_result/ -o ancestral_reconstruction/reconstruction_global_rpoB_pseudo -t new_tree/out_tree.nwk -m global -p {}".format(includePseudo)
        cmd14 = "./ancestral_reconstruction/reconstruction.py -i ancestral_reconstruction/new_result/ -o ancestral_reconstruction/reconstruction_global_operon_pseudo -t operon_tree/tree/out_tree.nwk -m global -p {}".format(includePseudo)
        print ("Executing command: {} \n {}\n".format(cmd13," reconstruction using rpoB"))
        os.system(cmd13)    
        print ("Executing command: {} \n {}\n".format(cmd14," reconstruction using operon tree"))
        os.system(cmd14) 
        ## show the tree
        cmd15 = "./ancestral_reconstruction/show_tree.py -g group.txt -i ancestral_reconstruction/reconstruction_global_operon_pseudo/gibberellin -a accession_to_common.txt -l a -o ancestral_reconstruction/Ryan_operon_pseudo"
        cmd16 = "./ancestral_reconstruction/show_tree.py -g group.txt -i ancestral_reconstruction/reconstruction_global_rpoB_pseudo/gibberellin -a accession_to_common.txt -l a -o ancestral_reconstruction/Ryan_rpoB_pseudo"
        print ("Executing command: {} \n {}\n".format(cmd15," showing operon tree"))
        os.system(cmd15)    
        print ("Executing command: {} \n {}\n".format(cmd16," showing rpoB tree"))
        os.system(cmd16) 
    else:
        cmd13 = "./ancestral_reconstruction/reconstruction.py -i ancestral_reconstruction/new_result/ -o ancestral_reconstruction/reconstruction_global_rpoB -t new_tree/out_tree.nwk -m global -p {}".format(includePseudo)
        cmd14 = "./ancestral_reconstruction/reconstruction.py -i ancestral_reconstruction/new_result/ -o ancestral_reconstruction/reconstruction_global_operon -t operon_tree/tree/out_tree.nwk -m global -p {}".format(includePseudo)
        print ("Executing command: {} \n {}\n".format(cmd13," reconstruction using rpoB"))
        os.system(cmd13)    
        print ("Executing command: {} \n {}\n".format(cmd14," reconstruction using operon tree"))
        os.system(cmd14) 
        ## show the tree
        cmd15 = "./ancestral_reconstruction/show_tree.py -g group.txt -i ancestral_reconstruction/reconstruction_global_operon/gibberellin -a accession_to_common.txt -l a -o ancestral_reconstruction/Ryan_operon"
        cmd16 = "./ancestral_reconstruction/show_tree.py -g group.txt -i ancestral_reconstruction/reconstruction_global_rpoB/gibberellin -a accession_to_common.txt -l a -o ancestral_reconstruction/Ryan_rpoB"
        print ("Executing command: {} \n {}\n".format(cmd15," showing operon tree"))
        os.system(cmd15)    
        print ("Executing command: {} \n {}\n".format(cmd16," showing rpoB tree"))
        os.system(cmd16) 

main()