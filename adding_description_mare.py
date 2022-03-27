#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 16:54:37 2022

@author: umut
"""
import sys
import csv
import pandas as pd
from os import listdir
from os.path import isfile, join

#print(sys.argv[0])

ddint = pd.read_csv('DrugBank_interactions.tsv', sep='\t')
df = pd.read_csv(sys.argv[1], '\t', names=["V1","V2","V3", "V4", "V5", "V6"])
df["V7"]='No'
df["V8"]='-'
df = df.reset_index()  # make sure indexes pair with number of rows
for index, row in df.iterrows():
    V3=row["V3"]
    V4=row["V4"]
    found=ddint.loc[(ddint['DrugBankID_1']==V3) & (ddint['DrugBankID_2']==V4)]
    if len(found.index)!=0:
        df.at[index,"V7"]='Yes'
        df.at[index,"V8"]=found['Description'].values[0]
          
       
    found=ddint.loc[(ddint['DrugBankID_1']==V4) & (ddint['DrugBankID_2']==V3)]
    if len(found.index)!=0:
        df.at[index,"V7"]='Yes'
        df.at[index,"V8"]=found['Description'].values[0]
         

print("Processed file : " +  sys.argv[1])           
df.to_csv(sys.argv[1], sep='\t', )        
            
        
        
            
            
        
        
        
    
    
    
 
