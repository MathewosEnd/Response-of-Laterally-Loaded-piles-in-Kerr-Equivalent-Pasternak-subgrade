from abaqus import *
from abaqusConstants import *
import visualization
import os
from numpy import savetxt, genfromtxt, array

#input Parameter 
# pData_file = 'pData.csv'
# pData = genfromtxt(pData_file, delimiter=',')
# eta_range = [0.001, 0.5, 1.0]

#Naming list
# baseName = 'FHFB_'
# instanceName_= ('SOIL-1', 'PILE-1')

#loop through results
workingDir = r'C:\Abaqus Thesis Files\Linearly Increasing Modulus\FxHFB\eta_zero'
os.chdir(workingDir)

pData_file = 'Parameters_LIModulus_FxHFB_eta_zero.csv'
pData = genfromtxt(pData_file, delimiter=',')

# Names that are used often
instanceName = ('SOIL-1', 'PILE-1')

for i in range(len(pData)):
    # Model name
    modelName = 'LIModulus_FxHFB_eta_zero_Model_' + str(i)
    jobName = modelName
    #
    varPile_u1 =[]
    node_coord = []
    # 
    #open Odb
    odb = session.openOdb(jobName + '.odb')
    # save coordinates of nodes on the pile center line
    regionPile_CL = odb.rootAssembly.instances[instanceName[1]].nodeSets['PFRONTEDGE_3']
    for s in range(len(regionPile_CL.nodes)):
        node_coord.append(regionPile_CL.nodes[s].coordinates[2])
    # obtain u1 from nodes on the pile center line and save
    stepLoading = odb.steps['Loading']
    lastFrame = stepLoading.frames[-1]
    var_u = lastFrame.fieldOutputs['U']
    varPile_u = var_u.getSubset(region=regionPile_CL)
    for k in range(len(varPile_u.values)):
        varPile_u1.append(varPile_u.values[k].data[0])
    output_matrix = array([node_coord, varPile_u1]).T
    savetxt(jobName + '.csv', output_matrix[output_matrix[:,0].argsort()], delimiter=',')
