import xml.etree.ElementTree as ET
import pandas as pd
import csv
import os 




#file name 
file_path = 'solution_test/1-rho_E-5-2-2.sol'

# Parse the XML file
tree = ET.parse(file_path)
root = tree.getroot()

#Extract values from tag  header information
header = root.find('header')
header_info = {
    'problemName': header.get('problemName'),
    'objectiveValue': header.get('objectiveValue'),
    'solutionTypeValue': header.get('solutionTypeValue'),
    'solutionTypeString': header.get('solutionTypeString'),
    'solutionStatusValue': header.get('solutionStatusValue'),
    'solutionStatusString': header.get('solutionStatusString'),
    'solutionMethodString': header.get('solutionMethodString'),
    'primalFeasible': header.get('primalFeasible'),
    'dualFeasible': header.get('dualFeasible'),
    'simplexIterations': header.get('simplexIterations'),
    'writeLevel': header.get('writeLevel')
}


#Extract values from tag quality metrics
quality = root.find('quality')
quality_metrics = {
    'epRHS': quality.get('epRHS'),
    'epOpt': quality.get('epOpt'),
    'maxPrimalInfeas': quality.get('maxPrimalInfeas'),
    'maxDualInfeas': quality.get('maxDualInfeas'),
    'maxPrimalResidual': quality.get('maxPrimalResidual'),
    'maxDualResidual': quality.get('maxDualResidual'),
    'maxX': quality.get('maxX'),
    'maxPi': quality.get('maxPi'),
    'maxSlack': quality.get('maxSlack'),
    'maxRedCost': quality.get('maxRedCost'),
    'kappa': quality.get('kappa')
}

#Extract values from tag linear constraints
constraints = []
for constraint in root.find('linearConstraints'):
    constraints.append({
        'name': constraint.get('name'),
        'index': constraint.get('index'),
        'slack': constraint.get('slack'),
    })

constraints_df = pd.DataFrame(constraints)



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

#convert type of columns
resultats = variables_df
resultats['value'] = resultats['value'].astype(float)
resultats['value'] = round(resultats['value']).astype(int)
resultats = resultats[resultats['value'] > 0]

# check 
print(resultats)


# # check
# print(variables_df['name'])

# file_csv = file_path[:-4]+'.csv'
# file_txt = file_path[:-4]+'.txt'


# variables_df.to_csv(file_csv,quoting=csv.QUOTE_NONNUMERIC,index=False,header=False)
# variables_df.to_csv(file_txt,index=False,header=False,sep='|',mode='a')


