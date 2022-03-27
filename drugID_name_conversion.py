import pandas
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

with open("drug_targeted.csv", 'r') as f:
    #for line in csv.reader(f):
        my_dict_2 = {rows[0]:rows[251] for rows in csv.reader(f)}  
        
for k, v in my_dict_2.items():
    line=v.split(';')
    newline=[]
    for j in line:
        temp=j.strip().lower()
        found=False
        for key, value in my_dict.items():
            if temp == value.strip().lower():
                newline.append(key)
                found=True
                break
        if not found:
            newline.append(j)
    newstr=""
    for part in newline:
        if part == newline[-1]:  
            newstr += part
        else:
            newstr += part + ';'
    my_dict_2[k]=newstr    
       
a_file = open("drug_targeted_conversion.csv", "w")

writer = csv.writer(a_file)
for key, value in my_dict_2.items():
    writer.writerow([key, value])

a_file.close()           
        
        
    
