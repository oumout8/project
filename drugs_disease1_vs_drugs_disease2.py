#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 12:09:01 2021

@author: umut
"""
import pandas as pd

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

#to find common drugs in PPI and CTD files (output pairs are from PPI file)
dictionary_ppi_intersection={}
for k, v in dictionary_ppi.items():
    line=v.split(';')
    newline=[]
    for j in line:
        found=False
        for key, value in dictionary_drug.items():
            if j == value:
                newline.append(j)
                found=True
                break     
        dictionary_ppi_intersection[k]=newline 
              
#filtering the PPI index with at least one value
filtered = {k: v for k, v in dictionary_ppi_intersection.items() if v}
dictionary_ppi_intersection.clear()
dictionary_ppi_intersection.update(filtered)   

disease2_input = str(input("Please Enter the disease name:"))

drug_disease2 = pd.read_csv('converted_filtered_drug_disease_after_conversion_with_disease_of_interest_without_dublicate.tsv', sep='\t')
disease_name2 = drug_disease2.loc[(drug_disease2['DiseaseName']==disease2_input)]
dictionary_drug2 = disease_name2.set_index('index_drug_disease')['ChemicalName'].to_dict()

drug_ppi2 = pd.read_csv('PPI_with_DrugIDs_converted_disease_name.tsv', sep='\t')
column_interested2=[col for col in drug_ppi2.columns if disease2_input in col]
selected_ppi2=drug_ppi2[drug_ppi2[disease2_input] == 1] 
dictionary_ppi2=selected_ppi2.set_index("index")["targeting.drugs.IDs"].to_dict()    
    
#to find common drugs in PPI and CTD files (output pairs are from CTD file)
dictionary_drug_intersection2={}
for k, v in dictionary_ppi2.items():
    line=v.split(';')
    for j in line:
        found=False
        for key, value in dictionary_drug2.items():
            if j == value:
                found=True
                break     
        dictionary_drug_intersection2[key]=value   

#to find common drugs in PPI and CTD files (output pairs are from PPI file)
dictionary_ppi_intersection2={}
for k, v in dictionary_ppi2.items():
    line=v.split(';')
    newline=[]
    for j in line:
        found=False
        for key, value in dictionary_drug2.items():
            if j == value:
                newline.append(j)
                found=True
                break     
        dictionary_ppi_intersection2[k]=newline 
              
#filtering the PPI index with at least one value
filtered = {k: v for k, v in dictionary_ppi_intersection2.items() if v}
dictionary_ppi_intersection2.clear()
dictionary_ppi_intersection2.update(filtered)

#To retrieve the common drugs in given disease pairs with disease 2 index
common_drugs_in_both_disease_2 = {}
for k, v in dictionary_drug_intersection.items():   
    found=False
    for key, value in dictionary_drug_intersection2.items():
        if v == value:
            found=True
            common_drugs_in_both_disease_2[key]=v

#To retrieve the common drugs in given disease pairs with disease 1 index
common_drugs_in_both_disease = {}
for k, v in dictionary_drug_intersection.items():   
    found=False
    for key, value in dictionary_drug_intersection2.items():
        if v == value:
            found=True
            common_drugs_in_both_disease[k]=v         

#To find different drugs between two diseases respect to disease 1
set1_disease1= set(dictionary_drug_intersection.items())
set2_disease1= set(common_drugs_in_both_disease.items())
not_common_drugs_in_both_disease_1= dict(set1_disease1-set2_disease1)

#To find different drugs between two diseases respect to disease 2
set1_disease2= set(dictionary_drug_intersection2.items())
set2_disease2= set(common_drugs_in_both_disease_2.items())
not_common_drugs_in_both_disease_2= dict(set1_disease2-set2_disease2)
          
'''          
#writing the result a file

with open('10.common_drugs_in_both_disease_with_index_disease_1.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['CTD FIle Index Disease 1','DrugIDs','Disease 1', 'Disease 2',])
    for key, val in common_drugs_in_both_disease.items():
        tsv_writer.writerow([key, val,disease_input ,disease2_input])


with open('11.common_drugs_in_both_disease_with_index_disease_2.tsv', 'wt') as out_file_2:
    tsv_writer_2 = csv.writer(out_file_2, delimiter='\t')
    tsv_writer_2.writerow(['CTD FIle Index Disease 2','DrugIDs','Disease 1', 'Disease 2',])
    for key, val in common_drugs_in_both_disease_2.items():
        tsv_writer_2.writerow([key, val,disease_input ,disease2_input])

with open('12.not_common_drugs_in_first_disease.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['CTD FIle Index Disease 1','DrugIDs','Disease 1', 'Disease 2',])
    for key, val in not_common_drugs_in_both_disease_1.items():
        tsv_writer.writerow([key, val,disease_input ,disease2_input])

with open('13.not_common_drugs_in_second_disease.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['CTD FIle Index Disease 1','DrugIDs','Disease 1', 'Disease 2',])
    for key, val in not_common_drugs_in_both_disease_2.items():
        tsv_writer.writerow([key, val,disease_input ,disease2_input])        
'''
import csv       
#appending the results to the same file
with open('10.common_drugs_in_both_disease_with_index_disease_1.tsv', "a") as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    for k, v in common_drugs_in_both_disease.items():
        tsv_writer.writerow([k,v,disease_input ,disease2_input])
out_file.close()        

with open('11.common_drugs_in_both_disease_with_index_disease_2.tsv', "a") as out_file_2:
    tsv_writer_2 = csv.writer(out_file_2, delimiter='\t')
    for k, v in common_drugs_in_both_disease_2.items():
        tsv_writer_2.writerow([k,v,disease_input ,disease2_input])
out_file.close()

with open('12.not_common_drugs_in_first_disease.tsv', "a") as out_file_2:
    tsv_writer_2 = csv.writer(out_file_2, delimiter='\t')
    for k, v in not_common_drugs_in_both_disease_1.items():
        tsv_writer_2.writerow([k,v,disease_input ,disease2_input])
out_file.close()
           
with open('13.not_common_drugs_in_second_disease.tsv', "a") as out_file_2:
    tsv_writer_2 = csv.writer(out_file_2, delimiter='\t')
    for k, v in not_common_drugs_in_both_disease_2.items():
        tsv_writer_2.writerow([k,v,disease_input ,disease2_input])
out_file.close()