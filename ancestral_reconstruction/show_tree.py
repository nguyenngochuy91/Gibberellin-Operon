#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Providing visualization of the Reconstruction. Grouping genomes
              into group color
    Start   : 05/08/2016
    End     : 05/08/2016
'''
from ete3 import *
import argparse
import os
from findParent_global import *

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Operon","-i", help="Operon file name")
    parser.add_argument("--Accession","-a", help="accession_to_common")
    parser.add_argument("--Group","-g", help="Color Grouping")
    parser.add_argument("--Image","-o", help="Output Image")
    parser.add_argument("--Locus","-l", help="gene mapping name (CYP114-h)")
    args = parser.parse_args()
    return args
    
# color group of genome
# create a dictionary, sadly this is done manually
def parse(file):
    color_dic={}
    group= open(file,'r')
    for line in group.readlines():
        key= line.split(':')[0]
        color=line.split(':')[1]
        value = color.split('\n')[0]
        color_dic[key]=value
    return color_dic

def accession_to_name(file):
    my_dic ={}
    infile = open(file,'r')
    for line in infile.readlines():
        line = line.replace('\r','')
        line = line.strip('\n')
        line = line.split(',')
        my_dic[line[0]] = line[1]
    return my_dic
if __name__ == "__main__":
    start = time.time()
    args = get_arguments()
    locus = args.Locus
    name_dic = accession_to_name(args.Accession)
    tree= Tree(args.Operon)
    mapping = args.Operon+'_mapping'
    infile = open(mapping,'r')
    dic={}
    for line in infile.readlines():
        line = line.strip()
        line = line.split('\t')
    for item in line:
        item = item.split(',')
        dic[item[1]]=item[0]
    color_list=['green','cyan','magenta','gray','yellow','orange',
               'red','lime','pink','blue','silver','maroon']
    gene_color_dic = {}
    for gene in dic:
        color = color_list.pop(0)
        gene_color_dic[gene]= color
    # file to write out about genome that has full othorlog CYP115
    outfile = open('full_CYP115.txt','w')
    # using the color dic to color group
    # color_dic = parse(args.Group)
    # fusion = open('potential_fusion.txt','w')
    pseudo = open('pseudo.txt','w')
    
    for node in tree.iter_descendants("postorder"):
        if not node.is_leaf():
            genes = list(node.initial)
            col = 1
            for gene in genes:
                if gene !="|":
                    if gene !="p":
                        gene_face = TextFace(gene)
                        gene_face.background.color = gene_color_dic[gene]           
                    else:
                        gene_face = TextFace(gene,fgcolor="white")
                        gene_face.background.color = "black"
                else:
                    gene_face = TextFace("  ")
                    gene_face.background.color = "white"
                node.add_face(gene_face,col,"branch-top")
                col+=1
                node.add_face(TextFace(" "),col,"branch-top")
                col+=1
            deletion_cost = (node.deletion).split('|')[1]
            dup_cost = (node.duplication).split('|')[1]
            split_cost = (node.split).split('|')[1]
            
            distances = [int(deletion_cost),int(dup_cost),int(split_cost)]
#            node.add_face(TextFace(node.initial), column=0, position = "branch-top")
            node.add_face(TextFace(distances), column=0, position = "branch-bottom")
            
            child1,child2 = node.get_children()
            color1 = child1.node_color
            color2 = child2.node_color
            if color1 == color2 and color1 !='mixed':
                node.add_features(node_color=color1)
            else:
                node.add_features(node_color='mixed')
        else:
            
                # print node.name, node.gene_block
            name = node.name.split("_")
            node.name =name_dic[node.name] 
            if locus in node.gene_block:
                outfile.write(node.name+'\n')
            color = 'red'
            node.add_features(node_color=color)
#            R = RectFace(1,2,color="RoyalBlue", label="123")
#            node.add_face(R)
            if "reference" in node.name:
                node.add_face(TextFace(node.name.replace("(reference)"," "),fgcolor = 'blue'), column =0, position ="aligned")
            elif "Erwinia_tracheiphila_PSU-1" in node.name:
                node.add_face(TextFace(node.name,fgcolor = 'gray'), column =0, position ="aligned")
            else:
                # detect those that dont have gene g (maybe because of fusion):
                if len(node.gene_block) ==0:
                    node.add_face(TextFace(node.name,fgcolor = 'red'), column =0, position ="aligned")
                else:
                    to_add =""
                        #fusion.write(node.name+'\n')
                    if 'p' in node.gene_block and locus not in node.gene_block:
                        to_add +='!'
                        pseudo.write(node.name+'\n')
                    elif 'p' not in node.gene_block and locus in node.gene_block:
                        to_add +='*'
                        pseudo.write(node.name+'\n')
                    if 'k' in node.gene_block:
                        to_add += '?'
                    node.add_face(TextFace(to_add+node.name), column =0, position ="aligned")
                    
            genes = list(node.gene_block)
            col = 1
            for gene in genes:
                if gene !="|":
                    if gene!="p":
                        gene_face = TextFace(gene)
                        gene_face.background.color = gene_color_dic[gene]      
                    else:
                        gene_face = TextFace(gene,fgcolor="white")
                        gene_face.background.color = "black"
                else:
                    gene_face = TextFace("  ")
                    gene_face.background.color = "white"
                node.add_face(gene_face,col,"aligned")
                col+=1
                node.add_face(TextFace(" "),col,"aligned")
                col+=1
            
            # node.dist = distance
        if node.node_color != 'mixed':
            nstyle = NodeStyle()
            nstyle["fgcolor"] = color
            # nstyle["vt_line_color"]=color
            # nstyle["hz_line_color"]=color
            node.set_style(nstyle)
    outfile.close()
    # fusion.close()
    pseudo.close()

    ### get the total cost for each event:
    # get the 2 children of the tree
    children= []
    # print "tree",tree
    for child in tree.get_children():
        children.append(child)
    # print "children",children
    deletion_cost1 = (children[0].deletion).split('|')[1]
    deletion_cost2 = (children[1].deletion).split('|')[1] 
    duplication_cost1 = (children[0].duplication).split('|')[1]
    duplication_cost2 = (children[1].duplication).split('|')[1]
    split_cost1 = (children[0].split).split('|')[1]
    split_cost2 = (children[1].split).split('|')[1]
    
    deletion_total = int(deletion_cost1) + int(deletion_cost2)
    duplication_total = int(duplication_cost1) + int(duplication_cost2)
    split_total = int(split_cost1)+int(split_cost2)
    # modify tree style for better visualization
    tree_style = TreeStyle()
    tree_style.show_leaf_name = False
    tree_style.min_leaf_separation = 5
    tree_style.extra_branch_line_type = 0
    tree_style.draw_guiding_lines=True
    tree_style.guiding_lines_type = 1
    cost= TextFace("Deletion count: "+str(deletion_total)+
                                '   '+"Duplication count: "+str(duplication_total)
                                +'   '+"Split count: "+ str(split_total),fsize =10,penwidth=2)
    cost.margin_top =5
    cost.margin_bottom = 5
    cost.margin_left = 5
    cost.margin_right = 5
    cost.background.color = 'LightGreen'
    tree_style.title.add_face(cost, column=1)

    dic= {'a': 'CYP115', 'b': 'CYP112', 'c': 'CYP114', 'd': 'Fd', 'e': 'SDR', 'f': 'CYP117', 'g': 'GGPS', 'h': 'CPS', 
    'i': 'KS', 'j': 'IDI','k':'GGPS2','p':'CYP115-pseudo/fragment'}
    mystring =''
    for item in sorted(dic):
        mystring += item+':'+dic[item]+'         '
    mystring = TextFace(mystring,fsize =10)
    mystring.margin_top =5
    mystring.margin_bottom = 5
    mystring.margin_left = 20
    mystring.margin_right = 20
    mystring.background.color = 'LightBlue'
    tree_style.title.add_face(mystring, column=2)
    
    # render the image
    tree.render(args.Image+'.png',dpi=1000,tree_style=tree_style)
    tree.render(args.Image+'.pdf',dpi=1000,tree_style=tree_style)
    # tree.render(args.Image+'.pdf',dpi=1000,tree_style=tree_style)
#    tree.show(tree_style=tree_style)



