import pprint, cPickle
import sys
import math
import numpy as np
import treeGDHeat
#sys.path.append("/home/asmariyaz/heatmap_codes")
import itertools
#import heatmap as HeatTree
import math
import get_probs_output

#def phylo_orderFunc(phylo_handle):
#    phylo_order=[i.strip("\n") for i in open(phylo_handle).readlines()]
#    return phylo_order

#def getTxtnames(txtname):
#    txtnames=[i.strip("\n") for i in txtname][1:]
#    return txtnames

# while traversing the event dict, insert 0 and nan and arrange all the scores in a list
# of list and store it as value in a dictionary with key as operon_event_name
# This list is later converted to a string and written to a csv file for error checking purposes
# and also used for plotting of the heat map in TreeGDHeat program.
def undoing(data,phylo_order):
    operon_dict = {}
    combine_dict ={}
    for operon, event in data.items():
        first_org = [tup[0] for tup in event.keys()]
        second_org = [tup[1] for tup in event.keys()]
        OrgInOperon= set(first_org + second_org)
        combination_whole = []
        splits_list_whole = []
        duplications_list_whole = []
        deletions_list_whole= []
        for org in phylo_order:
            splits_list_single =[org]
            duplications_list_single =[org]
            deletions_list_single = [org]
            combination_single = [org]
            for next_org in phylo_order:
                if (org in OrgInOperon) and (next_org in OrgInOperon) and (org != next_org):
                   if event.has_key((org, next_org)) or event.has_key((next_org, org)):
                      if event.has_key((org,next_org)):
                         key = event[(org, next_org)]
                      else:
                         key = event[(next_org, org)]
                      splitVal= key["splits"]  
                      duplicationVal = key["duplications"]
                      deletionsVal = key["deletions"]
                      if math.isnan(splitVal):
                         key["splits"] = 0
                      if math.isnan(duplicationVal):
                         key["duplications"] = 0
                      if math.isnan(deletionsVal):
                         key["deletions"] = 0
                      splits_list_single.append(key['splits'])
                      duplications_list_single.append(key['duplications'])
                      deletions_list_single.append(key['deletions'])
                      combination_single.append(key['splits']+key["duplications"]+key["deletions"])
                elif (org in OrgInOperon) and (next_org in OrgInOperon) and (org == next_org):
                     splits_list_single.append(0)
                     duplications_list_single.append(0)
                     deletions_list_single.append(0)
                     combination_single.append(0)
                elif (org not in OrgInOperon) or (next_org not in OrgInOperon):
                     splits_list_single.append(float("NaN"))
                     duplications_list_single.append(float("NaN"))
                     deletions_list_single.append(float("NaN"))

    ##        if operon == 'ybgIJKL-nei':
    ##           print deletions_list_single
            splits_list_whole.append(splits_list_single)
            duplications_list_whole.append(duplications_list_single)
            deletions_list_whole.append(deletions_list_single)
            if len(combination_single) >1:
                combination_whole.append(combination_single)
        operon_dict[operon +'_'+'splits']= splits_list_whole
        operon_dict[operon +'_'+'deletions']= deletions_list_whole
        operon_dict[operon +'_'+'duplications']= duplications_list_whole
        combine_dict[operon +'_'+'combination'] = combination_whole
#    print "operon_dict[0]",operon_dict[operon +'_'+'splits'][0]
#    print "combine_dict[0]",combine_dict[operon +'_'+'combination'][0]

    return operon_dict,combine_dict
    
def convert_str(phylo_order,listOfEvents,i):
    if i==0:
       phylo_order.insert(0,"Organism")
    header=",".join(phylo_order)+"\n"
    result=[header]

    for l in listOfEvents:
        combined = ','.join([str(i) for i in l])+ "\n"      
        result.append(combined)    
    return result
            
def writing_to_file(operon_event,result,csvDir):
    writing_filename=operon_event
    f=open(csvDir+ "/" +writing_filename + ".csv","w")
    f.writelines(result)
    return writing_filename
#######################################################################################

# Used to determine min and max of every operon's event's matrix
# These min and max values will be used as min and max of the colorbar in TreeGDHeat 
def minMax(L):
    result = list(itertools.chain(*L))
    #print "result",result
    min_val = min(filter(lambda t: not math.isnan(t), result))
    #print "min_val",min_val 
    max_val = max(filter(lambda t: not math.isnan(t), result))
    #print "max_val",max_val 
    return min_val, max_val

# Revised version for minMax, accomodate where all the thing is NaN for operon event
def minMax1(L):
    result = list(itertools.chain(*L))
    #print "result",result
    newResult = []
    for item in result:
        if item == item:
            newResult.append(item)
    #print "newResult",newResult
    if len(newResult) == 0: # all the thing in result was nan
        min_val = float('nan')
        max_val = float('nan')
    else:
        min_val = min(newResult)
        max_val = max(newResult)
    return min_val,max_val

#######################################################################################
def prep_event(event):
    result=[]
    for l in event:
        result.append(l[1:])
    return result


###              ###            ###
###Implementation for TreeHeat###
###              ###            ###
##pkl_file = open('event_dict_new.p.events.pik', 'rb')
##data1 = cPickle.load(pkl_file)
##phylo_order = phylo_orderFunc(phyloFile)
##txtnames = txtnames(txtnameFile)
##operon_dict=undoing(data1,phylo_order)
##i = 0
##full_max = 0
##full_min = 0
##combinedDir_path = HeatTree.determineDirectory()
##for operon,event in operon_dict.items():
##    result=convert_str(phylo_order, event,i)
##    full_len = prep_event(event)
##    min_num, max_num = minMax(full_len)
##    HeatTree.combineAll(max_num, min_num, full_len, txtnames, tree_path, operon, combinedDir_path)
##    f=writing_to_file(operon,result)
##    i+=1








#              ###            ###
#Implementation for TreeGDHeat###
#              ###            ###
def generateCombined(eventDict,legend_data,accession_order,organism_order,tree_path,csvDir,genomeDiagramsDir,combinedDir_path,tempDir):                
    eventZScoreFile,eventZScore = get_probs_output.run_main(eventDict, tempDir)
    operonEvents = cPickle.load(open(eventZScoreFile,'rb'))
    ##operonEvents = cPickle.load(eventDict)
    operon_dict,combine_dict = undoing(operonEvents,accession_order)
    i = 0
    full_max = 0
    full_min = 0
    gd_files_list = treeGDHeat.traverseAll(genomeDiagramsDir)
    #print "operon_dict",operon_dict
    for operon,event in operon_dict.items():
        result=convert_str(accession_order, event,i)
        full_len = prep_event(event)
        #print "event",event
        #print "operon", operon
        #print "full_len",full_len
        # min_num, max_num = minMax(full_len)
        min_num, max_num = minMax1(full_len)
        treeGDHeat.combineAll(max_num, min_num, full_len, organism_order, tree_path, operon, gd_files_list, combinedDir_path,
legend_data)
        f=writing_to_file(operon,result,csvDir)
        i+=1
    
    for operon,event in combine_dict.items():
        result=convert_str(accession_order, event,i)
        full_len = prep_event(event)
        #print "event",event
        #print "operon", operon
        #print "full_len",full_len
        # min_num, max_num = minMax(full_len)
        f=writing_to_file(operon,result,csvDir)
        i+=1
    

















'''
    eventZScoreFile,eventZScore = get_probs_output.run_main(eventDict, tempDir)
    operonEvents = cPickle.load(open(eventZScoreFile,'rb'))
    operon_dict=undoing(operonEvents,accession_order)
    i = 0
    full_max = 0
    full_min = 0
    gd_files_list = treeGDHeat.traverseAll(genomeDiagramsDir)

    for operon,event in operon_dict.items():
        result=convert_str(accession_order, event,i)
        full_len = prep_event(event)
        min_num, max_num = minMax(full_len)
        treeGDHeat.combineAll(max_num, min_num, full_len, organism_order, tree_path, operon, gd_files_list, combinedDir_path, legend_data)
        f=writing_to_file(operon,result,csvDir)
        i+=1
'''
