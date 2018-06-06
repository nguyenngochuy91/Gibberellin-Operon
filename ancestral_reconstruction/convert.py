#!/usr/bin/env python
''' Author  : Huy Nguyen, and David C.Ream
    Program : Given the operon directory, for each operons file, get the info about the gene in each genomes, map it
              into alphabet letter , get the gap, map gap to '|' and write to ouputfile.
    Start   : 05/04/2016
    End     : 05/05/2016
'''

import os
import argparse
import time
import uuid
import sys
# traverse and get the file
GPPS2 = "XOC_0086"
def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+f)
    return res

class readable_dir(argparse.Action):
    def __call__(self,parser, namespace, values, option_string=None):
        prospective_dir=values
        if not os.path.isdir(prospective_dir):
           try:
               os.mkdir(prospective_dir)
           except OSError:
               print (argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir)))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace,self.dest,prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--OperonDataDirectory","-i",action=readable_dir,help="This directory should contain files with gene name, start, stop, strand direction information for each genome.")
    parser.add_argument("--OutputDirectory","-o", help="Output of this program will be stored in the path supplied here. It will make a new directory if path given is valid or it will raise an error")
    parser.add_argument("--LocusName","-l",help="This contains info about the gene name (XOC_00..)")
    args = parser.parse_args()
    return args


def chk_output_directory_path(OutputDirectory,sessionID):
    if not os.path.exists(OutputDirectory + "_" + str(sessionID)):
        try:
           #os.mkdir(OutputDirectory + "_" + str(sessionID))
           return True
        except OSError:
           print ("Unable to create directory:", OutputDirectory)
           sys.exit()

# convert the file into dictionary with useful info    
def toDict(file,locus):
    infile = open(file,'r')
    map_code=''
    mapping =''
    dic_map ={}
    main_dic={}
    GPPS_file = open("GPPS2.txt","w")
    # create 3 main key
    first_line = infile.readline()
    map_code = first_line
    mapping = first_line.split('\t')[:-1] 
    for item in mapping:
        item_split = item.split(',')
        key = item_split[0]
        value = item_split[1]
        dic_map[key]=value # {'astA': 'a'}
    
    for line in infile.readlines():
        genome = line.split(':')[0] # (line.split(':')= ['NC_002696', '(astA,634700,635744,1)\t(astD,635730,637149,1)\t(astB,637145,638426,1)\t(astC,638435,639614,1)\t']
        main_dic[genome]={}
        # 3 distance sub dictionary
        main_dic[genome]['+1']={}
        main_dic[genome]['-1']={}
        main_dic[genome]['1']={}
        genes_string = line.split(':')[1]
        # to deal with each genes, they are in tuple, where first is the name of the gene, follow by the position, and the strand it is on
        # should consider 2 type of strand (so i can put a gap correctly
        genes_string = genes_string .split('\t')[:-1] # ['(astA,634700,635744,1)', '(astD,635730,637149,1)', '(astB,637145,638426,1)', '(astC,638435,639614,1)']
        genes_string = list(set(genes_string))
        if genome == "CP000133.1":
            print genes_string
        for item in genes_string:
            info= item.split(',') #['dppA', '402362', '400796', '+1','65.12']
            identity = float(info[4])
            locus_name = info[0]
            boolean = info[5]
            position=(int(info[1]),int(info[2]))
            position=(min(position),max(position))
            if locus_name == locus and  (boolean =="False" or identity<50):#dealing with pseudo gene/fragment of j
                main_dic[genome][info[3]][position]= 'p'
            else:
                if identity >=40 and locus_name != "XOC_0079":
                    main_dic[genome][info[3]][position]=dic_map[locus_name]
                if identity >=50 and locus_name == "XOC_0079":
                    main_dic[genome][info[3]][position]= 'g'
#                if identity <50 and locus_name == "XOC_0079": # pseudo gene of GGPS (d)
#                    main_dic[genome][info[3]][position]= 'o'
                if locus_name == "XOC_0086" and identity >=70:
                    main_dic[genome][info[3]][position]= 'k'
                GPPS_file.write(genome+":"+str(info)+"\n")
        if genome == "CP000133.1":
            print main_dic[genome]

    # print ("main_dic",main_dic)
            
    return (main_dic,map_code)

# from dic, create string
def toString(dic,map_code):
    wholestring=''
    wholestring+=map_code
    for genome in dic:
        string= genome + ':' # the string to be written
        substring = ''
        flag = False # check if it has a gene block
        for key in dic[genome]:
            if len(dic[genome][key])==0:
                continue
            else:
                myList=[]
                for position in dic[genome][key]:
                    myList.append(position)
                myList.sort()
                substring += dic[genome][key][myList[0]]
                for index in range(len(myList)-1):
                    dif = abs(myList[index+1][0] - myList[index][1]) 
                    if dif >800:
                        substring += '|'
                    else:
                        flag = True   
                    substring += dic[genome][key][myList[index+1]]
                substring += '|' # changing strand
        if flag:
            string += substring[:-1] # only add if there is a gene block
                
        string += '\n'
        wholestring += string
    return wholestring

        
if __name__ == "__main__":

    start = time.time()
    args = get_arguments()
    locus = args.LocusName
    sessionID = uuid.uuid1()
    condition = chk_output_directory_path(args.OutputDirectory,sessionID)
    if condition:
        outputsession = args.OutputDirectory
        try:
            os.mkdir(outputsession)
        except:
            print "dir already existed"
        res = traverseAll(args.OperonDataDirectory)
        for r in res:
            root,f = os.path.split(r)
            result= toDict(r,locus)
            wholestring = toString(result[0],result[1])
            outfile = open(outputsession+'/'+f,'w')
            outfile.write(wholestring)
            outfile.close()
    print (time.time() - start)