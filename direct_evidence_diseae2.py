#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 14:12:14 2021

@author: umut
"""

import csv

my_dict = {}

with open("11.common_drugs_in_both_disease_with_index_disease_2.tsv", 'r') as f:
    for line in f:
        items = line.split("\t")
        key, values = items[0], items[1]
        my_dict[key] = values


print (my_dict)

my_dict_2 = {}

with open("converted_filtered_drug_disease_after_conversion_with_disease_of_interest_without_dublicate.tsv", 'r') as f:
    my_dict_2 = {rows[0]:rows[6] for rows in csv.reader(f,delimiter = '\t', doublequote=False)} 
 
print (my_dict_2) 

for k, v in my_dict_2.items():   
    found=False
    for key, value in my_dict.items():
        if k == value:
            found=True
            my_dict[key]=v 

a_file = open("16.direct_evidence_disease_2.tsv", "w")
writer = csv.writer(a_file, delimiter="\t")
for key, value in my_dict.items():
    writer.writerow([key, value])

a_file.close()                      