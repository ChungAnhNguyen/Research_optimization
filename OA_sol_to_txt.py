import xml.etree.ElementTree as ET
import pandas as pd
import csv
import os 
import sys 

#file name 
filename = sys.argv[1]
file_path = './solution_test/'+filename

# Parse the XML file
tree = ET.parse(file_path)
root = tree.getroot()

#Extract values from tag variables
variables = []
text = []
for variable in root.find('variables'):
    text.append(variable.get('name')[1:])
    variables.append({
        'name': (variable.get('name'))[1:],
        # 'index': variable.get('index'),
        'value': variable.get('value'),
    })

variables_df = pd.DataFrame(variables).drop(index=0) #drop the first row that does not conaint important datas

#convert type of value's columns
resultats = variables_df
resultats['value'] = resultats['value'].astype(float)
resultats['value'] = round(resultats['value']).astype(int)
resultats = resultats[resultats['value'] > 0]
# check 
print(resultats)

#Convert to file 
file_csv = "./csv_file/"+filename[:-4]+'.csv'
file_txt = "./text_file/"+filename[:-4]+'.txt'
resultats.to_csv(file_csv,quoting=csv.QUOTE_NONNUMERIC,index=False,header=["Array","Frequent"])
resultats.to_csv(file_txt,index=False,header=False,sep='|',mode='a')