#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : find the representation of species
    Start   : 06/04/2017
    End     : 06/08/2017
'''


from ete3 import Tree
import argparse
def parser_code():

    parser = argparse.ArgumentParser(description='The purpose of this script to debias tree based on parameter')
    
    parser.add_argument("-f", "--full", help="full tree")
    
    parser.add_argument("-r", "--reduce", help="reduced tree")
                
    parser.add_argument("-o", "--output", help="output")
                                  
    return parser.parse_args()
def retrieveInfo(textFile):
    handle = open(textFile,"r")
    l = handle.readline()
    lines = handle.readlines()
    d = {}
    for line in lines:
        line = line.strip().split(":")
#        print (line)
        name = line[0]
        try:
            geneBlock = line[1]
        except:
            geneBlock = []
        d[name]   = []
        if "a" in geneBlock:
            d[name].append("a")
        if "k" in geneBlock:
            d[name].append("k") 
        d[name] = ",".join(d[name])
    return d
if __name__ == "__main__":
    args  = parser_code()
    infile = open("accession_to_common.txt","r")
    accessionName = {}
    infoDic       = retrieveInfo("ancestral_reconstruction/new_result/gibberellin")
#    print (infoDic)
    for line in infile.readlines():
        info = line.strip().split(",")
        accessionName[info[0]]= info[1]
    fullTree = Tree(args.full)
    partialTree = Tree(args.reduce)
    outfile = open(args.output,"w")
    representation = {}
    leavesFull = fullTree.get_leaves()
    allNames   = [leaf.name for leaf in partialTree]
    leavesPartial = partialTree.get_leaves()
    for leaf in fullTree:
        leafName = leaf.name
        if leafName not in allNames:
#            print ("leafName:",leafName)
            # only for those not show up, walk up until hit one that is shown, check its sibbling, or sibling children
            ancestors = leaf.get_ancestors()
            for ancestor in ancestors:
                subtreeLeaf = ancestor.get_leaves()
                # check if the leaf appear in all Names
                found = False
                for subLeaf in subtreeLeaf:
                    name = subLeaf.name
                    if name in allNames:
#                        print ("representation:",name)
                        found = True
                        # add this name to dictionary
                        try:
                            fromInfoDic = infoDic[name]
                        
                        except:
                            fromInfoDic = ''
                        try:
                            leafInfoDic = infoDic[leafName]
                        except:
                            leafInfoDic = ''
                        if (name,accessionName[name],fromInfoDic) in representation:
                            representation[(name,accessionName[name],fromInfoDic)].append((leafName,accessionName[leafName],leafInfoDic))
                        else:
                            representation[(name,accessionName[name],fromInfoDic)] = [(leafName,accessionName[leafName],leafInfoDic)]
                        break
                if found:
                    break
    # write into outfile
    swapFile = open("swap_{}".format(args.output),"w")
    for name in representation:
        string = "{}:".format(name)
        nameList = ["{}".format(item) for item in representation[name]]
        string+=",".join(nameList)      
        string+="\n"
        outfile.write(string)
        if name[2]=="":
            check = False
            for n in representation[name]:
                if n[2] =="k" or n[2]=="a":
                    swapFile.write("{}:{}\n".format(name,n))
    outfile.close()
    swapFile.close()