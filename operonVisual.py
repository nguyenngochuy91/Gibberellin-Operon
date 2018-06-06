#!/usr/bin/env python

import os
import sys
import reportlab.lib
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation,SeqFeature
import ntpath
from matplotlib.colors import colorConverter 
from matplotlib.colors import rgb2hex
import argparse
import pickleToCSV
import uuid
import time
from homolog4 import *

# written by David, modified by Huy
# TODO: parallelize this code, it takes forever to execute and single threaded on an
# embarrassingly parallel problem

## Traverses the genome information directory
def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+f)
    return res


#def reading_phyloFile(file_handle):
#    list_lines=[i.strip("\n") for i in file_handle]
#    return list_lines

def reading_optFile(file_handle):
    list_lines=[i.strip("\n") for i in open(file_handle,'rU').readlines()]
    return list_lines

##
## parses the genome information file
##def listToDict(list_lines):
##    result_dict={}
##    for l in list_lines:
##        if l.startswith("NC_"):
##            list_tab=l.split("\t")
##            accession=list_tab[0]
##            geneName=list_tab[4]
##            start=list_tab[10]
##            end=list_tab[11]
##            strand_direction=list_tab[12]
##            if result_dict.has_key(accession)==False:
##               result_dict[accession]=[(geneName,start,end,strand_direction)]
##            else:
##               dummy_dict={}
##               list1=result_dict[accession]
##               list2=[(geneName,start,end,strand_direction)]
##               list1.extend(list2)
##               dummy_dict[accession]=list1
##               result_dict.update(dummy_dict)         
##    return result_dict
##

def gdVisualizationDict(infile):
    result = {}
    for line in [i.strip() for i in open(infile).readlines()]:
            threshold = 300 
            hlog = Homolog.from_blast(line)
            accession = hlog.accession()
            start = hlog.start()
            end = hlog.stop()
            strand = hlog.strand()
            # gene_name = hlog.blast_annotation()
            locus = hlog.query_locus()
            identity = hlog.percent_ident()
            boolean = abs(hlog.align_subject_start()-hlog.align_subject_stop())<threshold
            # print "identity",identity
            if accession in result.keys():
                result[accession].append((locus, start, end, strand,identity,boolean))
            else:
                result.update({accession:[(locus, start, end, strand,identity,boolean)]})
    # print "result_gd", result
    return result




#
def allReverseStranded(strand_list):
    condition = all(x == strand_list[0] for x in strand_list)
    if condition == True:
       if strand_list[0] == str(-1) or strand_list[0] == -1:
           condition = True
       else:
           condition = False
    return condition

def changeAllStrandDirection(geneStartStopStrandDirection):
    ## find '-1' convert it to '+1'
    for i,tup in enumerate(geneStartStopStrandDirection):
        new_tup = (tup[0],tup[2],tup[1],'+1',tup[4],tup[5])
        geneStartStopStrandDirection[i] = new_tup
    return geneStartStopStrandDirection

def SplittingList(geneStartStopStrandDirection, accession_num,i,j):
    splitList = []
    for ele in geneStartStopStrandDirection:
           
        if not splitList:
            splitList.append([ele]) # final is empty so add the first element in a list
        # check the last item of the last sublist added to final and compare to our current element
        elif abs(int(splitList[-1][-1][i]) - int(ele[j])) > 500 or splitList[-1][-1][3] != ele[3]:
        # if it does not meet the requirement, add it to final in a new list
        #     if accession_num == 'NC_003047' and 'atp' in geneStartStopStrandDirection[0][0]:
        #        print accession_num+":"+splitList[-1][-1][i] + "-" + ele[j] + "diff:" + str(abs(int(splitList[-1][-1][1]) - int(ele[2])))
        #     if accession_num == 'NC_000913' and 'atp' in geneStartStopStrandDirection[0][0]:
        #        print accession_num+":"+ splitList[-1][-1][i] + "-" + ele[j] + "diff:" + str(abs(int(splitList[-1][-1][1]) - int(ele[2])))   
               
             splitList.append([ele])
        else:
            # else add it to the last sublist
            splitList[-1].append(ele)   
    return splitList

def findMajority(geneStartStopStrandDirection , acc):
    i = 2
    j = 1
    splitList = SplittingList(geneStartStopStrandDirection , acc, i ,j)
    maxLen = 0
    maxList = []
    for subList in splitList:
        if (len(subList) > maxLen):
            maxLen = len(subList)
            maxList = subList        
    strand = maxList[0][3]
    if strand == '-1' or strand == -1:
        strandedness = 'negative'
    else:
        strandedness = 'positive'
    return strandedness, splitList, maxList

def findEqual(geneStartStopStrandDirection):    
    if (len(geneStartStopStrandDirection) % 2) == 0:
        positiveNum = 0
        negativeNum = 0
        for tup in geneStartStopStrandDirection:
            if tup[3] == '-1' or tup[3] == -1:
                negativeNum+=1
            else:
                positiveNum+=1
        if positiveNum == negativeNum:
            condition = True
        else:
            condition = False
    else:
         condition = False
    return condition    

def changeStrandDirection(eleTup):
    newEleTup = (eleTup[0],eleTup[2],eleTup[1],'+1',eleTup[4],eleTup[5])
    return newEleTup

def changeStrandDirectionOption(eleTup,direction):
    newEleTup = (eleTup[0],eleTup[2],eleTup[1],direction,eleTup[4],eleTup[5]) 
    return newEleTup

def handle_strandedness(parsed_Dict):
    strandCorrDict = {}
    for accession_num, geneStartStopStrandDirection in parsed_Dict.items():
        if geneStartStopStrandDirection:
            strand_list = []
            ## check if all the genes are on the reverse strand
            for tup in geneStartStopStrandDirection:
                strand_list.append(tup[3])
            condition = allReverseStranded(strand_list)
            if condition == True:
               Changed_geneStartStopStrandDirection = changeAllStrandDirection(geneStartStopStrandDirection)
         ##      if accession_num == 'NC_000913' and 'atp' in geneStartStopStrandDirection[0][0] :
         ##         print Changed_geneStartStopStrandDirection
               i = 1
               j = 2
               Split_geneStartStopStrandDirection = SplittingList(Changed_geneStartStopStrandDirection,accession_num,i,j)
               combinedList = []
               for l in Split_geneStartStopStrandDirection:
                   for ele in reversed(l):
                       combinedList.append(ele)
               strandCorrDict[accession_num] = combinedList
            else:
               ## check if equal number of genes point in both directions
               equal = findEqual(geneStartStopStrandDirection) 
               if equal == False:
                   direction, SplitList, maxList = findMajority(geneStartStopStrandDirection,accession_num) 
                   if direction == 'negative':
          ##            if accession_num == 'NC_003047' and 'atp' in SplitList[0][0][0]:
          ##               print "SplitList:",SplitList
          ##               print "maxList:",maxList
                      reversedMajorityList = [] 
                      for l in SplitList:
                          if l == maxList:
                             for ele in reversed(l):
                                 new_ele = changeStrandDirection(ele)
                                 reversedMajorityList.append(new_ele)
                          else:
                             for ele in reversed(l):
                                 strandDirection = ele[3]
                                 if strandDirection == str(1) or strandDirection == '+1' or strandDirection == 1:
                                    direction = '-1'
                                 else:
                                    direction = '+1'
                                 new_ele = changeStrandDirectionOption(ele,direction)       
                                 reversedMajorityList.append(new_ele)
           ##           if accession_num == 'NC_003047' and 'atp' in SplitList[0][0][0]:
           ##              print 'reversedMajorityList:',reversedMajorityList
                      
                      strandCorrDict[accession_num] = reversedMajorityList      
                   else:
                      strandCorrDict[accession_num] = geneStartStopStrandDirection
               else:
                  strandCorrDict[accession_num] = geneStartStopStrandDirection  
    #print "strandCorrDict",strandCorrDict
    return strandCorrDict

#assigns a reportlab color to each gene
def geneToColor(geneList):
    #print "geneList",geneList
    idToColorDict_reportLab={}
    idToColorDict_matplotlib={}

    Colorlist1 = [colors.red,colors.blue,colors.orange,colors.plum,colors.green,colors.gray,colors.black,colors.magenta,
                                     colors.chartreuse,colors.sandybrown,colors.darkturquoise,colors.indigo,colors.chocolate,
                                     colors.aquamarine]
    Colorlist2=["red","blue",
                   "orange","plum",
                   "green","gray",
                   "black","magenta",
                   "chartreuse","sandybrown",
                   "darkturquoise","indigo",
                   "chocolate","aquamarine"]
    for i,g in enumerate(geneList):
        #print "i",i
        #print "g",g
        idToColorDict_reportLab[g]=Colorlist1[i]
        idToColorDict_matplotlib[g]=Colorlist2[i]
    return idToColorDict_reportLab,idToColorDict_matplotlib

# Traverses the parsed file(stored as a Dict) and draws a genome diagram
def drawGenomeDiag(result_Dict,phylo_order,split_distance,filename,OutputDirectory):
    
    #print "result_dic", result_Dict
    operon_name = filename.split('.')[0]
    # print "operon_name",operon_name

    gd_diagram =GenomeDiagram.Diagram(filename)
    geneList=[]
    maxNoGenes = 0
    for k,v_l in result_Dict.items():
        for v in v_l:
            geneList.append(v[0])
        if maxNoGenes < len(geneList):
           maxNoGenes = len(geneList)
    set_gene = set(geneList)
    # mapping gene to a alphabet letter for my algorithm
    sort_gene = list(set_gene)
    sort_gene.sort()

    idToColorDict_reportLab,idToColorDict_matplotlib=geneToColor(set_gene)
    for rec in phylo_order:

        cnt = 0
        gd_track_for_features= gd_diagram.new_track(1,name=rec,greytrack=False,start=0,end=200)
        gd_feature_set=gd_track_for_features.new_set()
        if result_Dict.has_key(rec):
           start=0
           maxPosList =[]
           minPosList = []
           tmp_list = result_Dict[rec]
           # print "tmp_list",tmp_list
           maxStrandList = []
           ##Code for changing the strand
           '''forw_count = 0
           revs_count = 0
           revrs_flag = False
           for l in tmp_list:
               strandDirection=l[3]
               if strandDirection == str(1) or strandDirection == '+1' or strandDirection == 1:
                   forw_count = forw_count + 1
               else:
                   revs_count = revs_count + 1
           tmp_list1 = []
           if(forw_count < revs_count):
               for l in reversed(tmp_list):
                   strandDirection=l[3]
                   temp = []
                   temp.append(l[0])
                   temp.append(l[1])
                   temp.append(l[2])
                   if strandDirection == str(1) or strandDirection == '+1' or strandDirection == 1:
                       temp.append('-1')
                   elif strandDirection == str(-1) or strandDirection == -1:
                       temp.append('+1')
                   else:
                       temp.append(l[3])
                   tmp_list1.append(temp)
           else:
               tmp_list1 = tmp_list'''
               
           for l in tmp_list:

               geneName=l[0]
               startPos=int(l[1])
               endPos=int(l[2])
               strandDirection=l[3]
               identity = l[4]
               boolean  = l[5]
               # adding the percent identity to show pseudo/ fragment
               string = l[0] + ','+ str(startPos)+ ','+ str(endPos) + ','+ str(strandDirection)+','+str(identity)+','+str(boolean)+'\t'

               if start==0:
                  start,maxPosList,maxStrandList,minPosList = findFeature(start,gd_feature_set,idToColorDict_reportLab,geneName,maxPosList,
                                                               minPosList,maxStrandList,endPos,startPos,strandDirection)
                  cnt+=1
               else:
                  if strandDirection != maxStrandList[-1]:
                     start = start + 10
                     start,maxPosList,maxStrandList,minPosList = findFeature(start,gd_feature_set,idToColorDict_reportLab,
                                                                  geneName,maxPosList,minPosList,maxStrandList,endPos,startPos,strandDirection)
                     cnt+=2
                  else:
                     diff = startPos - maxPosList[-1]
                     if abs(diff) < split_distance or (maxPosList[-1] - endPos == 0 and minPosList[-1] - startPos == 0):                     
                        start,maxPosList,maxStrandList,minPosList = findFeature(start,gd_feature_set,idToColorDict_reportLab,
                                                                     geneName,maxPosList,minPosList,maxStrandList,endPos,startPos,strandDirection)
                        cnt+=1
                     else:
                         start = start + 10
                         start,maxPosList,maxStrandList,minPosList = findFeature(start,gd_feature_set,idToColorDict_reportLab,
                                                                      geneName,maxPosList,minPosList,maxStrandList,endPos,startPos,strandDirection)
                         cnt+=2

    gd_diagram.draw(format="linear", pagesize="A4",fragments=1,start=0,end=180,track_size=0.7)
    gd_diagram.write(OutputDirectory + "/" + filename + ".png","PNG")
    return idToColorDict_matplotlib
    
def findFeature(start,gd_feature_set,idToColorDict_reportLab,geneName,maxPosList,minPosList,maxStrandList,endPos,startPos,strandDirection):
    if strandDirection == str(1) or strandDirection == '+1' or strandDirection == 1:
       fz= SeqFeature(FeatureLocation(start,end=start+10,strand=+1))
    elif strandDirection == str(-1) or strandDirection == -1:
       fz= SeqFeature(FeatureLocation(start,end=start+10,strand=-1))
    else:
       print type(strandDirection)
    gd_feature_set.add_feature(fz,sigil="BIGARROW",color = idToColorDict_reportLab[geneName],
                               label=False,name=geneName,label_position="start")                
    start=start+10
    maxPosList.append(endPos)
    minPosList.append(startPos)
    maxStrandList.append(strandDirection)
    return start,maxPosList,maxStrandList,minPosList

class readable_dir(argparse.Action):
    def __call__(self,parser, namespace, values, option_string=None):
        prospective_dir=values
        if not os.path.isdir(prospective_dir):
           try:
               os.mkdir(prospective_dir)
           except OSError:
               print argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace,self.dest,prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

# reading in a csv file with Accession ids and respective organism names
def reading_MappingFile(mappingFile):
    reader = open(mappingFile,'r')
    accession = []
    organism = []
    for line in reader.readlines():
        line = line.strip().split(',')
        accession.append(line[0])
        organism.append(line[1])
    return accession, organism

#def removeTrailingBackSlashes(dirpath):
#    if dirpath.endswith("/"):
#       cleanedpath = dirpath[:-1]
#    else:
#       cleanedpath = dirpath
#    return cleanedpath

def makeSubfolder(output_directory_path, subfolder_name, sessionID):
    path = output_directory_path + "/" + subfolder_name + "_" + str(sessionID)
    os.makedirs(path)
    return path

def makeSubfolder2(output_directory_path, subfolder_name):
    path = output_directory_path + "/" + subfolder_name
    os.makedirs(path)
    return path


def chk_output_directory_path(OutputDirectory,sessionID):
    if not os.path.exists(OutputDirectory + "_" + str(sessionID)):
        try:
           #os.mkdir(OutputDirectory + "_" + str(sessionID))
           return True
        except OSError:
           print "Unable to create directory:", OutputDirectory
           sys.exit()
             

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--OperonDataDirectory","-n",action=readable_dir,help="This directory should contain files with gene name, start, stop, strand direction information.")
    parser.add_argument("--MappingFile","-m",help="This file must contain Genbank Accession numbers and its mapping to organism names in order of the newick tree nodes, csv format required.")
    parser.add_argument("--NewickTree","-t",help="This file must be a newick formatted tree.")
    parser.add_argument("--EventsDict","-e",type=file, help="This file contains all the operons' events and pairwise distances.")
    parser.add_argument("--splitDistance","-d", type=int ,default=500,help="Splitting distance")
    parser.add_argument("--OutputDirectory","-o", help="Output of this program will be stored in the path supplied here. It will make a new directory if path given is valid or it will raise an error")
    args = parser.parse_args()
    return args

if __name__ == "__main__":

    start = time.time()

    all_colors = list(reportlab.lib.colors.getAllNamedColors().items())
    args = get_arguments()
    sessionID = uuid.uuid1()
    condition = chk_output_directory_path(args.OutputDirectory,sessionID)
    if condition:
       outputsession = args.OutputDirectory # + "_" + str(sessionID)    
       '''    
       OutputGenomeDiagDirectory = makeSubfolder(outputsession,"genome-diagrams",sessionID)
       OutputCSVDirectory = makeSubfolder(outputsession,"operon-event-matrices",sessionID)
       OutputTreeGDHeatDirectory = makeSubfolder(outputsession,"tree-gd-heat-diagrams",sessionID)
       TempDirectory = makeSubfolder(outputsession, "temporary-files",sessionID)
       '''
       OutputGenomeDiagDirectory = makeSubfolder2(outputsession,"genome-diagrams")
       OutputCSVDirectory = makeSubfolder2(outputsession,"operon-event-matrices")
       OutputTreeGDHeatDirectory = makeSubfolder2(outputsession,"tree-gd-heat-diagrams")
       TempDirectory = makeSubfolder2(outputsession, "temporary-files")
       mapping = args.MappingFile
       print mapping
       accession_order,organism_order = reading_MappingFile(mapping)
       res = traverseAll(args.OperonDataDirectory)
       split_distance = args.splitDistance
       legendData = {}
       #print "*****"
       print "Results can be found in the following directory:", outputsession
       #print "*****"

       for r in res:
           root,f = os.path.split(r)
           ##listLines = reading_optFile(r)
           ##result_dict = listToDict(listLines)
           # print f;
           result_dict = gdVisualizationDict(r)
           changedStrandedness = handle_strandedness(result_dict)
           # if "caiTABCDE" in r: 
           #    print "caiTABCDE", r
           #    print changedStrandedness
           #print "result_dict", changedStrandedness
           idToColorDict_matplotlib = drawGenomeDiag(changedStrandedness,accession_order,split_distance,f,OutputGenomeDiagDirectory)
           ##print "legend_data",ntpath.basename((r.split("/")[4]).split(".")[0])
           
           operonName =  ntpath.basename(r)
           #print operonName;
           legendData[operonName.split(".")[0]] = idToColorDict_matplotlib
       pickleToCSV.generateCombined(args.EventsDict,legendData,accession_order,organism_order,args.NewickTree,OutputCSVDirectory,OutputGenomeDiagDirectory,OutputTreeGDHeatDirectory,TempDirectory)   
       
       
    print time.time() - start
