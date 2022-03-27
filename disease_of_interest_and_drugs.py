#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 25 01:19:03 2021

@author: umut
"""

import pandas as pd
import csv

disease_input = str(input("Please Enter the disease name:"))

drug_disease = pd.read_csv('converted_filtered_drug_disease_after_conversion_with_disease_of_interest_without_dublicate.tsv', sep='\t')
disease_name = drug_disease.loc[(drug_disease['DiseaseName']==disease_input)]
dictionary_drug = disease_name.set_index('index_drug_disease')['ChemicalName'].to_dict()

drug_ppi = pd.read_csv('PPI_with_DrugIDs_converted_disease_name.tsv', sep='\t')
column_interested=[col for col in drug_ppi.columns if disease_input in col]
selected_ppi=drug_ppi[drug_ppi[disease_input] == 1] 
dictionary_ppi=selected_ppi.set_index("index")["targeting.drugs.IDs"].to_dict()    

#to find common drugs in PPI and CTD files (output pairs are from CTD file)
dictionary_drug_intersection={}
for k, v in dictionary_ppi.items():
    line=v.split(';')
    for j in line:
        found=False
        for key, value in dictionary_drug.items():
            if j == value:
                found=True
                break     
        dictionary_drug_intersection[key]=value   
  
with open('18.drugs_with_CTD_index_disease_of_interest.tsv', 'a') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['CTD FIle Index','DrugIDs','Disease ',])
    for key, val in dictionary_drug_intersection.items():
        tsv_writer.writerow([key, val,disease_input])        