#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 15:11:03 2021

@author: umut
"""

import csv
import itertools

my_dict = {}

with open("DrugBank_names.tsv", 'r') as f:
    for line in f:
        items = line.split("\t")
        key, values = items[0], items[1]
        my_dict[key] = values


print (my_dict)

my_dict_2 = {}

with open("filtered_drug_disease.tsv",  "rt") as f:
        my_dict_2 = {rows[0]:rows[1] for rows in csv.reader(f,delimiter = '\t', doublequote=False)}  
        
#my_dict_2=dict(itertools.islice(my_dict_2.items(), 1000))


for k, v in my_dict_2.items():   
    found=False
    for key, value in my_dict.items():
        if v.strip().lower() == value.strip().lower():
            found=True
            my_dict_2[k]=key 
    if not found:
        my_dict_2[k]=v        
      
 
my_dict_3 = {}

with open("DrugBank_names.tsv", 'r') as f:
    for line in f:
        items = line.split("\t")
        key, values = items[0], items[2]
        my_dict_3[key] = values


#print (my_dict_3)

for k, v in my_dict_3.items():
    line=v.split(';')
    for j in line:
        temp=j.strip().lower()
        found=False
        for key, value in my_dict_2.items():
            if temp == value.strip().lower():
                found=True
                my_dict_2[key]=k
        if not found:
            my_dict_2[key]=value    

        
a_file = open("drug_disease_conversion.tsv", "w")
writer = csv.writer(a_file, delimiter="\t")
for key, value in my_dict_2.items():
    writer.writerow([key, value])

a_file.close()          


      