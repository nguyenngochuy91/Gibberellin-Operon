#!/usr/bin/python

# Copyright(C) 2015 Ashish Jain
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

regulon_url = 'http://regulondb.ccg.unam.mx/menu/download/datasets/files/OperonSet.txt'
regulon_outfolder = 'regulonDB/'
regulon_default_gene_block_file = regulon_outfolder + 'gene_block_names_and_genes.txt'
regulon_experimental_only = True
tree_marker_gene = "rpob"
tree_outfolder = "tree/"
BLAST_database_folder = 'db/'
BLAST_format_protein = 'True'
gene_block_query_outfile = 'gene_block_query.fa'
refrence_organism = 'reference'
BLAST_outfolder = 'blast_result/'
BLAST_parse_outfolder = 'blast_parse/'
filter_BLAST_parse_outfolder = 'optimized_gene_blocks/'
event_distance_outfolder = 'gene_block_distance_matrices/'
visualization_outfolder = "visualization/"
newick_tree = 'out_tree.nwk'
accession_to_common = 'accession_to_common.csv'
phylo_order_new = 'phylo_order_new.txt'
event_types = ['deletions','splits','duplications']
