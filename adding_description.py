#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 16:54:37 2022

@author: umut
"""
import csv
import pandas as pd
from os import listdir
from os.path import isfile, join

ddint = pd.read_csv('DrugBank_interactions.tsv', sep='\t')

mypath= "/home/umut/Desktop/Thesis_Project/Oumout/Files/Separation_Score_Intersect"
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

for i, f in enumerate(onlyfiles):
    df = pd.read_csv(mypath+'/'+ f, '\t')#, names=["Drug1","Drug2","d_A", "d_B", "s_AB"])
    df["Interaction"]='No'
    df["Description"]='-'
    df = df.reset_index()  # make sure indexes pair with number of rows
    for index, row in df.iterrows():
        V3=row["Drug1"]
        V4=row["Drug2"]
        found=ddint.loc[(ddint['DrugBankID_1']==V3) & (ddint['DrugBankID_2']==V4)]
        if len(found.index)!=0:
            df.at[index,"Interaction"]='Yes'
            df.at[index,"Description"]=found['Description'].values[0]
            #print("I have found: " + row["Drug1"] +" " +row["Drug2"] +" " + found['Description'].values[0])
       
        found=ddint.loc[(ddint['DrugBankID_1']==V4) & (ddint['DrugBankID_2']==V3)]
        if len(found.index)!=0:
            df.at[index,"Interaction"]='Yes'
            df.at[index,"Description"]=found['Description'].values[0]
            #print("I have found: " + row["Drug1"] +" " +row["Drug2"] +" " + found['Description'].values[0])

    print(str(i)+" "+"Processed file : " + f )           
    df.to_csv(mypath+'/'+ f, sep='\t', )        
            
        
        
            
            
        
        
        
    
    
    
 
