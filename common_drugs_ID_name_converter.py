#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 15:34:44 2021

@author: umut
"""

import csv

my_dict = {}

with open("10.common_drugs_in_both_disease_with_index_disease_1.tsv", 'r') as f:
    for line in f:
        items = line.split("\t")
        key, values = items[0], items[2]
        my_dict[key] = values


print (my_dict)

my_dict_2 = {}

with open("DrugBank_names.tsv", 'r') as f:
    for line in f:
        items = line.split("\t")
        key, values = items[0], items[1]
        my_dict_2[key] = values


print (my_dict_2)

for k, v in my_dict_2.items():   
    found=False
    for key, value in my_dict.items():
        if k == value:
            found=True
            my_dict[key]=v 

a_file = open("14.names_of_common_drugs_in_both_disease_with_index_disease_1.tsv", "w")
writer = csv.writer(a_file, delimiter="\t")
for key, value in my_dict.items():
    writer.writerow([key, value])

a_file.close()          