#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 13:18:16 2022

@author: umut
"""

import pandas as pd
import csv

df = pd.read_csv("DrugBank_targets.tsv", sep='\t')
df.drop(columns=['TargetGeneSymbol','Action'], inplace=True)
mydict=df.groupby('DrugBankID')['TargetUniProtID'].apply(list).to_dict()

a_file = open("20.DrugBank_targets_concatenated.tsv", "w")
writer = csv.writer(a_file, delimiter="\t")
for key, value in mydict.items():
    writer.writerow([key, value])

a_file.close() 


my_dict = {}

with open("19.drugs_with_CTD_index_disease_of_interest_with_index.tsv", 'r') as f:
    for line in f:
        items = line.split("\t")
        key, values = items[0], items[2]
        my_dict[key] = values


print (my_dict)
my_dict = {key.replace('"', ''):val.replace('"', '' ) for key, val in my_dict.items()}

my_dict_2 = {}

with open("20.DrugBank_targets_concatenated.tsv", 'r') as f:
    for line in f:
        items = line.strip().split("\t")
        k, v = items[0], items[1]
        my_dict_2[k] = v


print (my_dict_2)

for key, values in my_dict.items():   
    found=False
    for k, v in my_dict_2.items():
        if k == values:
            found=True
            my_dict[key]=v 

a_file = open("21.DrugBank_targets_concatenated_after_conversion.tsv", "w")
writer = csv.writer(a_file, delimiter="\t")
for key, value in my_dict.items():
    writer.writerow([key, value])

a_file.close()             