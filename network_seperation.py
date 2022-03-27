#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 11:01:43 2022

@author: umut
"""

import networkx as nx
import numpy as np
import csv
import optparse
import sys
import pandas as pd
import os

def read_network(network_file):
    """
    Reads a network from an external file.

    * The edgelist must be provided as a tab-separated table. The
    first two columns of the table will be interpreted as an
    interaction gene1 <==> gene2

    * Lines that start with '#' will be ignored
    """

    G = nx.Graph()
    for line in open(network_file,'r'):
        # lines starting with '#' will be ignored
        if line[0]=='#':
            continue
        # The first two columns in the line will be interpreted as an
        # interaction gene1 <=> gene2
        line_data   = line.replace('"', '').strip().split('\t')
        node1 = line_data[2]
        node2 = line_data[3]
        G.add_edge(node1,node2)

    print ("\n> done loading network:")
    print ("> network contains %s nodes and %s links" %(G.number_of_nodes(),
                                                       G.number_of_edges()))
    
    return G

# =============================================================================
def read_gene_list(gene_file):
    """
    Reads a list genes from an external file.

    * The genes must be provided as a table. If the table has more
    than one column, they must be tab-separated. The first column will
    be used only.

    * Lines that start with '#' will be ignored
    """

    genes_set = set()
    for line in open(gene_file,'r'):
        # lines starting with '#' will be ignored
        if line[0]=='#':
            continue
        # the first column in the line will be interpreted as a seed
        # gene:
        line_data = line.strip().split('\t')
        gene      = line_data[0]
        genes_set.add(gene)

    print ("\n> done reading genes:")
    print ("> %s genes found in %s" %(len(genes_set),gene_file))

    return genes_set

# =============================================================================
def remove_self_links(G):

    sl = G.selfloop_edges()
    G.remove_edges_from(sl)



# =============================================================================
def get_pathlengths_for_single_set(G,given_gene_set):
    
    """
    calculate the shortest paths of a given set of genes in a
    given network. The results are stored in a dictionary of
    dictionaries:
    all_path_lenghts[gene1][gene2] = l
    with gene1 < gene2, so each pair is stored only once!

    PARAMETERS:
    -----------
        - G: network
        - gene_set: gene set for which paths should be computed

    RETURNS:
    --------
        - all_path_lenghts[gene1][gene2] = l for all pairs of genes
          with gene1 < gene2

    """ 

    # remove all nodes that are not in the network
    all_genes_in_network = set(G.nodes())
    gene_set = given_gene_set & all_genes_in_network

    all_path_lenghts = {}
    
    # calculate the distance of all possible pairs
    for gene1 in gene_set:
        if gene1 not in all_path_lenghts:
            all_path_lenghts[gene1] = {}
        for gene2 in gene_set:
            if gene1 < gene2:
                try:
                    l = nx.shortest_path_length(G, source=gene1, target=gene2)
                    all_path_lenghts[gene1][gene2] = l
                except:
                    continue

    return all_path_lenghts



# =============================================================================
def get_pathlengths_for_two_sets(G,given_gene_set1,given_gene_set2):
    
    """
    calculate the shortest paths between two given set of genes in a
    given network. The results are stored in a dictionary of
    dictionaries: all_path_lenghts[gene1][gene2] = l with gene1 <
    gene2, so each pair is stored only once!

    PARAMETERS:
    -----------
        - G: network
        - gene_set1/2: gene sets for which paths should be computed

    RETURNS:
    --------
        - all_path_lenghts[gene1][gene2] = l for all pairs of genes
          with gene1 < gene2

    """ 

    # remove all nodes that are not in the network
    all_genes_in_network = set(G.nodes())
    gene_set1 = given_gene_set1 & all_genes_in_network
    gene_set2 = given_gene_set2 & all_genes_in_network

    all_path_lenghts = {}
    
    # calculate the distance of all possible pairs
    for gene1 in gene_set1:
        if gene1 not in all_path_lenghts:
            all_path_lenghts[gene1] = {}
        for gene2 in gene_set2:
            if gene1 != gene2:
                try:
                    l = nx.shortest_path_length(G, source=gene1, target=gene2)
                    if gene1 < gene2:
                        all_path_lenghts[gene1][gene2] = l
                    else:
                        if gene2 not in all_path_lenghts:
                            all_path_lenghts[gene2] = {}
                        all_path_lenghts[gene2][gene1] = l
                except:
                    continue

    return all_path_lenghts


# =============================================================================
def calc_single_set_distance(G,given_gene_set):

    """
    Calculates the mean shortest distance for a set of genes on a
    given network    
    

    PARAMETERS:
    -----------
        - G: network
        - gene_set: gene set for which distance will be computed 

    RETURNS:
    --------
         - mean shortest distance 

    """
    


    # remove all nodes that are not in the network, just to be safe
    all_genes_in_network = set(G.nodes())
    gene_set = given_gene_set & all_genes_in_network

    # get the network distances for all gene pairs:
    all_path_lenghts = get_pathlengths_for_single_set(G,gene_set)

    all_distances = []

    # going through all gene pairs
    for geneA in gene_set:

        all_distances_A = []
        for geneB in gene_set:

            # I have to check which gene is 'smaller' in order to know
            # where to look up the distance of that pair in the
            # all_path_lengths dict
            if geneA < geneB:
                if geneB in all_path_lenghts[geneA]:
                    all_distances_A.append(all_path_lenghts[geneA][geneB])
            else:
                if geneA in all_path_lenghts[geneB]:
                    all_distances_A.append(all_path_lenghts[geneB][geneA])

        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)

    # calculate mean shortest distance
    mean_shortest_distance = np.mean(all_distances)

    return mean_shortest_distance


# =============================================================================
def calc_set_pair_distances(G,given_gene_set1,given_gene_set2):

    """
    Calculates the mean shortest distance between two sets of genes on
    a given network
    
    PARAMETERS:
    -----------
        - G: network
        - gene_set1/2: gene sets for which distance will be computed 

    RETURNS:
    --------
         - mean shortest distance 

    """

    # remove all nodes that are not in the network
    all_genes_in_network = set(G.nodes())
    gene_set1 = given_gene_set1 & all_genes_in_network
    gene_set2 = given_gene_set2 & all_genes_in_network

    # get the network distances for all gene pairs:
    all_path_lenghts = get_pathlengths_for_two_sets(G,gene_set1,gene_set2)

    all_distances = []

    # going through all pairs starting from set 1 
    for geneA in gene_set1:

        all_distances_A = []
        for geneB in gene_set2:

            # the genes are the same, so their distance is 0
            if geneA == geneB:
                all_distances_A.append(0)
                
            # I have to check which gene is 'smaller' in order to know
            # where to look up the distance of that pair in the
            # all_path_lengths dict
            else:
                if geneA < geneB:
                    try:
                        all_distances_A.append(all_path_lenghts[geneA][geneB])
                    except:
                        pass

                else:
                    try:
                        all_distances_A.append(all_path_lenghts[geneB][geneA])
                    except:
                        pass


        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)

    # going through all pairs starting from disease B
    for geneA in gene_set2:

        all_distances_A = []
        for geneB in gene_set1:

            # the genes are the same, so their distance is 0
            if geneA == geneB:
                all_distances_A.append(0)

            # I have to check which gene is 'smaller' in order to know
            # where to look up the distance of that pair in the
            # all_path_lengths dict
            else:
                if geneA < geneB:
                    try:
                        all_distances_A.append(all_path_lenghts[geneA][geneB])
                    except:
                        pass
                        
                else:
                    try:
                        all_distances_A.append(all_path_lenghts[geneB][geneA])
                    except:
                        pass

        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)


    # calculate mean shortest distance
    mean_shortest_distance = np.mean(all_distances)

    return mean_shortest_distance

# LOADING NETWORK and DRUGS GENES
    #
    # --------------------------------------------------------
   
listofdiseases=["Asthma","Lung Diseases, Obstructive","Carcinoma, Non-Small-Cell Lung","Schizophrenia","Alzheimer Disease","Arthritis, Rheumatoid","Diabetes Mellitus, Type 2","Obesity","Hypertension", "Cardiovascular Diseases", "Atherosclerosis", "Bone Marrow Neoplasms", "Leukemia, Lymphocytic, Chronic, B-Cell","Lymphoma","Liver Neoplasms","Neuroblastoma", "Breast Neoplasms","Prostatic Neoplasms","Leukemia, Myeloid, Acute", "Urinary Bladder Neoplasms", "Glioblastoma","Adenocarcinoma","Colonic Neoplasms","Colorectal Neoplasms","Melanoma","Stomach Neoplasms","Pancreatic Neoplasms","Carcinoma, Squamous Cell"]    
drugs = pd.read_csv('22.drugs_with_CTD_index_disease_of_interest_with_protein_target.tsv', sep='\t')
ppin=  pd.read_csv('PPI_with_DrugIDs_converted_disease_name.tsv', sep='\t')
    
for z in range(len(listofdiseases)):
    ## Is just for the test, comment the following line
    #z<-1
    print("Starting Disease1: Number " +str(z)+" - " +listofdiseases[z])
    for y in range(z+1,len(listofdiseases)):
        ## Is just for the test, comment the following line
        #y<-z+1
        disease1=listofdiseases[z]
        disease2=listofdiseases[y]
        # print(paste(disease1, disease2,sep="-"))    
    
        #disease1="Asthma"
        #disease2="Carcinoma, Non-Small-Cell Lung"
    
        
        ppin_d1d2= ppin.loc[(ppin[disease1]==1) & (ppin[disease2]==1)]
    
        ppin_d1d2.to_csv('temporal_file.tsv', sep='\t', ) 

        # read network
        G  = read_network("temporal_file.tsv")  
        # get all genes
        all_genes_in_network = set(G.nodes())
    
        list_drug1=[]
        list_drug2=[]
        list_d_A=[]
        list_d_B=[]
        list_d_AB=[]
        list_s_AB=[]
        my_df=drugs[drugs["Disease "]==disease1][["DrugIDs","UniProtID"]]
        my_df2=drugs[drugs["Disease "]==disease2][["DrugIDs","UniProtID"]]
        my_df = my_df.reset_index() 
        my_df2= my_df2.reset_index()# make sure indexes pair with number of rows
        for index, row in my_df.iterrows():
            drug=row[1]
            proteins=row[2]
            genec=list(proteins.split(", "))
            genes_A = set(genec) & all_genes_in_network
            d_A = calc_single_set_distance(G,genes_A)
            d_A= np.nan_to_num(d_A)
            for index, row in my_df2.iterrows():
                drug2=row[1]
                proteins2=row[2]
                gened=list(proteins2.split(", "))
                genes_B = set(gened) & all_genes_in_network
                d_B = calc_single_set_distance(G,genes_B)
                d_B = np.nan_to_num(d_B)
                d_AB = calc_set_pair_distances(G,genes_A,genes_B)
                d_AB = np.nan_to_num(d_AB)
                s_AB = d_AB - (d_A + d_B)/2.
                list_drug1.append(drug)
                list_drug2.append(drug2)
                list_d_A.append(d_A)
                list_d_B.append(d_B)
                list_d_AB.append(d_AB)
                list_s_AB.append(s_AB)
        df = pd.DataFrame(list(zip(list_drug1, list_drug2, list_d_A, list_d_B,list_d_AB,list_s_AB, )), columns =['Drug1', 'Drug2','d_A','d_B','d_AB','s_AB'])
        os.chdir('/home/umut/Desktop/Thesis_Project/Oumout/Files/Separation_Score_Intersect')         
        df.to_csv(disease1+"_"+disease2, index=False, sep='\t')
        os.chdir('/home/umut/Desktop/Thesis_Project/Oumout/Files') 
    

    
    









    

