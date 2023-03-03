"""
Abaqus automation script for an elastic circular pile embeded in elastic soil under lateral loading.

This script attomates modeling of pile and elastic media with linearly increasing
elastic modulus with Es(0) = 0.001 (eta aprox. equal to 0). 
"""
#----------------------------------------------------------------------------------------------------------------
#                                    import neccesary packages
#----------------------------------------------------------------------------------------------------------------
from abaqus import *
from abaqusConstants import *
import mesh
import os
from numpy import pi, sin, cos, genfromtxt, Linspace
#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#                                    FUNCTION used to build and run model
#----------------------------------------------------------------------------------------------------------------
def createPart_Soil(modelName, partName, domainDia, domainDepth):
    #create the soil part
    partName = partName[0]
    s = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE) 
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(domainDia/2.0, 0.0))
    s.Line(point1=(-domainDia/2.0, 0.0), point2=(domainDia/2.0, 0.0))
    s.autoTrimCurve(curve1=g[2], point1=(0.0, domainDia/2.0))
    p = mdb.models[modelName].Part(name=partName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models[modelName].parts[partName]
    p.BaseSolidExtrude(sketch=s, depth=domainDepth)
    s.unsetPrimaryObject()
    p = mdb.models[modelName].parts[partName]
    del mdb.models[modelName].sketches['__profile__']

def createPilecut_Soil(modelName, partName, pileDia, pileLength):
    #create a cut in the soil for the pile
    partName = partName[0]
    p = mdb.models[modelName].parts[partName]
    f, e = p.faces, p.edges
    t = p.MakeSketchTransform(sketchPlane=f[3], sketchUpEdge=e[2], sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0))
    s1 = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=100, gridSpacing=2.0, transform=t)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p = mdb.models[modelName].parts[partName]
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, pileDia/2.0))
    s1.Line(point1=(0.0, pileDia/2.0), point2=(0.0, -pileDia/2.0))
    s1.autoTrimCurve(curve1=g[4], point1=(pileDia/2.0, 0.0))
    p.CutExtrude(sketchPlane=f[3], sketchUpEdge=e[2], sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s1, depth=pileLength, flipExtrudeDirection=OFF)
    s1.unsetPrimaryObject()
    del mdb.models[modelName].sketches['__profile__']

def createPart_Pile(modelName, partName, pileDia, pileLength):
    #create the soil part
    partName = partName[1]
    s = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE) 
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(pileDia/2.0, 0.0))
    s.Line(point1=(-pileDia/2.0, 0.0), point2=(pileDia/2.0, 0.0))
    s.autoTrimCurve(curve1=g[2], point1=(0.0, pileDia/2.0))
    p = mdb.models[modelName].Part(name=partName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models[modelName].parts[partName]
    p.BaseSolidExtrude(sketch=s, depth=pileLength)
    s.unsetPrimaryObject()
    p = mdb.models[modelName].parts[partName]
    del mdb.models[modelName].sketches['__profile__']

def createDatumplane(modelName, partName, datumPlane, datumDepth):
    p = mdb.models[modelName].parts[partName]
    datum = p.DatumPlaneByPrincipalPlane(principalPlane=datumPlane, offset=datumDepth)
    idDatum = datum.id
    return idDatum

def createPartion(modelName, partName, idDatum):
    p = mdb.models[modelName].parts[partName]
    c = p.cells[:]
    d = p.datums
    p.PartitionCellByDatumPlane(datumPlane=d[idDatum], cells=c)

def ExtrudPartition(modelName, partName):
    p = mdb.models[modelName].parts[partName]
    c = p.cells[:]
    e, d1 = p.edges, p.datums
    pickedEdges =(e[14], e[23])
    p.PartitionCellByExtrudeEdge(line=e[6], cells=c, edges=pickedEdges, sense=FORWARD)

def CreateRefPoint(modelName, setname, x, y, z):
    a = mdb.models[modelName].rootAssembly
    refPoint = a.ReferencePoint (point=(x, y, z))
    r = a.referencePoints
    refPoint_Loc = r.findAt((x, y, z), )
    refPoints1=(refPoint_Loc, )
    a.Set(referencePoints=refPoints1, name=setname)
    return refPoint, refPoint_Loc

def createPSet_cellWhole(modelName, partName, setName_cellWhole):
    #create a set for cells
    #Note No partitions are considered
    p = mdb.models[modelName].parts[partName]
    c = p.cells[:] #only one cell exists
    p.Set(cells=c, name=setName_cellWhole)

def createPset_cell(modelName, partName, setName, x ,y ,z):
    #create a set for cells
    #considers partitions
    cell = ()
    p = mdb.models[modelName].parts[partName]
    c = p.cells 
    cellLoc = c.findAt((x, y, z), )
    cell = cell + (c[cellLoc.index:cellLoc.index+1], )
    p.Set(cells=cell, name=setName)

def createIset_cell(modelName, instanceName, setName, x ,y ,z):
    #create a set for cells
    #considers partitions
    cell = ()
    a = mdb.models[modelName].rootAssembly
    c = a.instances[instanceName].cells 
    cellLoc = c.findAt((x, y, z), )
    cell = cell + (c[cellLoc.index:cellLoc.index+1], )
    iRegion = a.Set(cells=cell, name=setName)
    return iRegion

def createPSet_face(modelName, partName, setName, x, y, z):
    #create a set for faces
    face = ()
    p = mdb.models[modelName].parts[partName]
    f = p.faces
    faceLoc = f.findAt((x, y, z), )
    face = face + (f[faceLoc.index:faceLoc.index+1], )
    p.Set(faces=face, name=setName)

def createISet_face(modelName, instanceName, setName, x, y, z):
    #create a set for faces
    face = ()
    a = mdb.models[modelName].rootAssembly
    f = a.instances[instanceName].faces
    faceLoc = f.findAt((x, y, z), )
    face = face + (f[faceLoc.index:faceLoc.index+1], )
    a.Set(faces=face, name=setName)

def createPSet_edge(modelName, partName, setName, x, y, z):
    #create a set that includes the pile part
    edge = ()
    p = mdb.models[modelName].parts[partName]
    e = p.edges
    edgeLoc = e.findAt((x, y, z), )
    edge = edge + (e[edgeLoc.index:edgeLoc.index+1], )
    p.Set(edges=edge, name=setName)

def createISet_edge(modelName, instanceName, setName, x, y, z):
    #create a set that includes the pile part
    edge = ()
    a = mdb.models[modelName].rootAssembly
    e = a.instances[instanceName].edges
    edgeLoc = e.findAt((x, y, z), )
    edge = edge + (e[edgeLoc.index:edgeLoc.index+1], )
    a.Set(edges=edge, name=setName)

def createPSet_vertex(modelName, partName, setName, x ,y ,z):
    vertex = ()
    p = mdb.models[modelName].parts[partName]
    v = p.vertices
    vertexLoc = v.findAt((x, y, z), )
    vertex = vertex + (v[vertexLoc.index:vertexLoc.index+1], )
    p.Set(vertices=vertex, name=setName)

def createISurface_region(modelName, instanceName, surfName, x, y, z):
    #create a set for faces
    face = ()
    a = mdb.models[modelName].rootAssembly
    s = a.instances[instanceName].faces
    faceLoc = s.findAt((x, y, z), )
    face = face + (s[faceLoc.index:faceLoc.index+1], )
    surfRegion = a.Surface(side1Faces=face, name=surfName)
    return surfRegion

def defineProperties_soil(modelName, materialName, sectionName, Es0, EsL, nu_s, pileLength):
    #define properties for soil
    materialName = materialName[0]
    sectionName = sectionName[0]
    #material properties
    mdb.models[modelName].Material(name=materialName)
    mdb.models[modelName].materials[materialName].Elastic(temperatureDependency=ON, table=((Es0, nu_s, 0.0), (EsL, nu_s, pileLength)))
    #section property
    mdb.models[modelName].HomogeneousSolidSection(name=sectionName, material=materialName, thickness=None)
    #section property
    mdb.models[modelName].HomogeneousSolidSection(name=sectionName, material=materialName, thickness=None)


def assignProperties_soil(modelName, partName, sectionName, setName_cellWhole):
    #assgin properties
    sectionName = sectionName[0]
    p = mdb.models[modelName].parts[partName]
    region = p.sets[setName_cellWhole]
    p.SectionAssignment(region=region, sectionName=sectionName, offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

def defineProperties_pile(modelName, materialName, sectionName, Ep, nu_p):
    #define properties for pile
    materialName = materialName[1]
    sectionName = sectionName[1]
    #material properties
    mdb.models[modelName].Material(name=materialName)
    mdb.models[modelName].materials[materialName].Elastic(table=((Ep, nu_p), ))
    #section property
    mdb.models[modelName].HomogeneousSolidSection(name=sectionName, material=materialName, thickness=None)

def assignProperties_pile(modelName, partName, sectionName, setName_cellWhole):
    #assgin properties
    sectionName = sectionName[1]
    p = mdb.models[modelName].parts[partName]
    region = p.sets[setName_cellWhole]
    p.SectionAssignment(region=region, sectionName=sectionName, offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

def createAssembly(modelName, partName, instanceName):
    #create assembly
    for i in range(len(partName)):
        a = mdb.models[modelName].rootAssembly
        a.DatumCsysByDefault(CARTESIAN)
        p = mdb.models[modelName].parts[partName[i]]
        a.Instance(name=instanceName[i], part=p, dependent=OFF)

def createStep(modelName):
    #create step
    #name is fixed
    mdb.models[modelName].StaticStep(name='Loading', previous='Initial')

def createOutputReq(modelName, fieldoutputReqName):
    #create field output request
    mdb.models[modelName].FieldOutputRequest(name=fieldoutputReqName, createStepName='Loading', variables=('U', 'S'), frequency=LAST_INCREMENT)
    del mdb.models[modelName].fieldOutputRequests['F-Output-1']

def createConstraint_tie(modelName, mainSurface, secondSurface, conName):
    # create a tie constraint between pile and soil
    mdb.models[modelName].Tie(name=conName, main=mainSurface, secondary=secondSurface, positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

def createConstraint_Coupling(modelName, controlSet, secondarySurface, constraintName):
    a = mdb.models[modelName].rootAssembly
    region1=a.sets[controlSet]
    region2=a.surfaces[secondarySurface]
    mdb.models[modelName].Coupling(name=constraintName, controlPoint=region1, surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

def assignBC_BaseFx(modelName, instanceName, setName, fxName):
    # fix the base of the domain
    a = mdb.models[modelName].rootAssembly
    region = a.instances[instanceName].sets[setName]
    mdb.models[modelName].DisplacementBC(name=fxName, createStepName='Initial', region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

def assignBC_lateralFx(modelName, instanceName, setName, fxName):
    # fix the lateral extent of the domain
    a = mdb.models[modelName].rootAssembly
    region = a.instances[instanceName].sets[setName]
    mdb.models[modelName].DisplacementBC(name=fxName, createStepName='Initial', region=region, u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)    

def assignBC_Ysymm(modelName, instanceName, setName, fxName):
    # apply symmetric bc to the XZ-plane
    a = mdb.models[modelName].rootAssembly
    region = a.instances[instanceName].sets[setName]
    mdb.models[modelName].YsymmBC(name=fxName, createStepName='Initial', region=region, localCsys=None)

def assignPileHeadLatLoad(modelName, instanceName, surfName, loadName, loadMag):
    # assign lateral load to the pile head
    a = mdb.models[modelName].rootAssembly
    e1 = a.instances[instanceName].edges
    v1 = a.instances[instanceName].vertices
    region = a.surfaces[surfName]
    mdb.models[modelName].SurfaceTraction(name=loadName, createStepName='Loading', region=region, magnitude=loadMag, directionVector=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), distributionType=UNIFORM, field='', localCsys=None)

def assignPileHeadMoment(modelName, setName, momentName, momentMag):
    a = mdb.models[modelName].rootAssembly
    region = a.sets[setName]
    mdb.models[modelName].Moment(name=momentName, createStepName='Loading', region=region, cm2=momentMag, distributionType=UNIFORM, field='', localCsys=None)

def assingDTempField(modelName, instanceName, setName):
    mdb.models[modelName].ExpressionField(name='AnalyticalField-1', localCsys=None, description='', expression=' Z ')
    a = mdb.models[modelName].rootAssembly
    region = a.instances[instanceName].sets[setName]
    mdb.models[modelName].Temperature(name='DummyTemp', createStepName='Initial', region=region, distributionType=FIELD, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, field='AnalyticalField-1', magnitudes=(1.0, ))

def assignMeshControlWEDGE(modelName, instanceName, x, y, z):
    # assign wedge meshe geometery to the pile and soil below the pile
    a = mdb.models[modelName].rootAssembly
    c = a.instances[instanceName].cells 
    region = c.findAt((x, y, z), )
    a.setMeshControls(regions=(region, ), elemShape=WEDGE, technique=SWEEP)

def createSeed_EdgeNum(modelName, instanceName, seedNumber, x, y, z):
    # seed the circular edge of the soil and pile by number
    a = mdb.models[modelName].rootAssembly
    e = a.instances[instanceName].edges
    edge = e.findAt((x, y, z), ) 
    a.seedEdgeByNumber(edges=(edge, ), number=seedNumber, constraint=FINER)

def createSeed_EdgeBias_2(modelName, instanceName, seedSizeMin, seedSizeMax, x, y, z):
    #bias seed direction 1
    # TODO: use a bias method that changes direction based on adge orentiation
    a = mdb.models[modelName].rootAssembly
    e = a.instances[instanceName].edges
    edge = e.findAt((x, y, z), ) 
    a.seedEdgeByBias(biasMethod=SINGLE, end2Edges=(edge, ), minSize=seedSizeMin, maxSize=seedSizeMax, constraint=FINER)

def createSeed_EdgeBias_1(modelName, instanceName, seedSizeMin, seedSizeMax, x, y, z):
    # bias seed direction 2
    a = mdb.models[modelName].rootAssembly
    e = a.instances[instanceName].edges
    edge = e.findAt((x, y, z), ) 
    a.seedEdgeByBias(biasMethod=SINGLE, end1Edges=(edge, ), minSize=seedSizeMin, maxSize=seedSizeMax, constraint=FINER)

def createSeed_global(modelName, instanceName, seedSizeGlob):
    # global seed 
    # used for the vertical edges of the pile and surrounding soil
    a = mdb.models[modelName].rootAssembly
    partInstances =(a.instances[instanceName[0]], a.instances[instanceName[1]], )
    a.seedPartInstance(regions=partInstances, size=seedSizeGlob, deviationFactor=0.1, minSizeFactor=0.1)

def createMesh(modelName, instanceName):
    #create mesh
    a = mdb.models[modelName].rootAssembly
    partInstances =(a.instances[instanceName[0]], a.instances[instanceName[1]], )
    a.generateMesh(regions=partInstances)

def assignElemType_Hex(modelName, instanceName, x, y ,z):
    # assign hex element types
    # second order elements with reduced integration 
    # hybrid formulation used when Poission's ratio exceeds or equal to 0.45
    elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=C3D20RH, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
    elemType4 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
    a = mdb.models[modelName].rootAssembly
    c = a.instances[instanceName].cells 
    region = c.findAt((x, y, z), )
    if nu_s < 0.45:
        a.setElementType(regions=(region, ), elemTypes=(elemType1, elemType3, elemType4))
    else:
        a.setElementType(regions=(region, ), elemTypes=(elemType2, elemType3, elemType4))

def assignElemType_Wedge(modelName, instanceName, x, y ,z):
    # assign wedge element types
    # second order elements with reduced integration 
    # hybrid formulation used when Poission's ratio exceeds or equal to 0.45
    elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=C3D15H, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
    elemType4 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
    a = mdb.models[modelName].rootAssembly
    c = a.instances[instanceName].cells 
    region = c.findAt((x, y, z), )
    if nu_s < 0.45:
        a.setElementType(regions=(region, ), elemTypes=(elemType1, elemType3, elemType4))
    else:
        a.setElementType(regions=(region, ), elemTypes=(elemType1, elemType2, elemType4))

def createJob(jobName, modelName):
    # create job
    mdb.Job(name=jobName, model=modelName, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT, numCpus=6, numDomains=6, numGPUs=2)

def submitJob(jobName):
    # submmit job
    mdb.jobs[jobName].submit(consistencyChecking=OFF)
    mdb.jobs[jobName].waitForCompletion()

#----------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
#                                      Model Building
#----------------------------------------------------------------------------------------------------------------

#input parameters
#set Working directory
workingDir = r'C:\Abaqus Thesis Files\Linearly Increasing Modulus\FHFB\eta_zero'
os.chdir(workingDir)

pData_file = 'Parameters_LIModulus_FHFB_eta_zero.csv'
pData = genfromtxt(pData_file, delimiter=',')

# Names that are used often
partName = ('partSoil', 'partPile')
sectionName = ('sectSoil', 'sectPile')
materialName = ('elasticSoil', 'elasticPile')
instanceName = ('SOIL-1', 'PILE-1')
fieldoutputReqName = 'F-Output-2'

for i in range(len(pData)):
    # Model name
    modelName = 'LIModulus_FHFB_eta_zero_Model_' + str(i)
    # Parameters
    pileRad = pData[i, 0]
    slendernessRatio = pData[i, 1]
    #
    eta = pData[i, 2]
    nu_s = pData[i, 3]
    Ep_mr = pData[i, 4]
    nu_p = pData[i, 5]
    Ep = pData[i, 6]
    #
    pileHeadLoad = pData[i, 7]
    pileHeadMoment = pData[i, 8] 
    #
    # Computed Parameters
    #
    # Material Property Related
    mr = Ep/Ep_mr
    EsL = (mr*slendernessRatio)/(1.0-eta)
    if eta < 0.001: # for when eta is set as 0.0, the inquality is used to account for error that might occure from equating floating points
        Es0 = 0.001
    else:
        Es0 = eta*EsL
    #
    pileDia = 2.0*pileRad
    pileLength = pileRad * slendernessRatio
    #
    # Domian size related
    # domainDepth = max{pileLength+15.0m, z_factor*pileLength}
    z_factor = 2.0
    domainDepth = z_factor*pileLength
    if domainDepth-pileLength < 15.0:
        domainDepth = pileLength + 15.0
    #
    # domain radius = max{2.0*pileLength, 30.0m} 
    domainDia = 4.0*pileLength
    if domainDia < 60.0:
        domainDia = 60.0
    #
    # Mesh related
    elemSize_outer = 5.0
    seedSizeMax_R = elemSize_outer
    seedSizeMin_R = 0.1
    #
    if pileLength > 15.0:
        seedSizeGlob = 0.75
    elif pileLength <= 15.0 and pileLength >= 10.0:
        seedSizeGlob = 0.5
    else:
        seedSizeGlob = 0.25
    #
    seedSizeMin_Z = seedSizeGlob
    seedSizeMax_Z = 2.5*seedSizeGlob
    seedSizeMax_Z_1 = seedSizeGlob
    seedSizeMin_Z_1 = 0.2
    #
    perimeter = 0.5*pi*domainDia
    seedNumber = int(round(perimeter/(elemSize_outer*2.0)))
    #
    # pile head load
    pileHeadTraction = pileHeadLoad/((pi*pileDia**2)/4.0)
    #
    #
    #create model
    mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT)
    jobName = modelName
    #create soil geometry
    createPart_Soil(modelName, partName, domainDia, domainDepth)
    createPilecut_Soil(modelName, partName, pileDia, pileLength)
    idDatum = createDatumplane(modelName, partName[0], datumPlane = XYPLANE, datumDepth = pileLength)
    createPartion(modelName, partName[0], idDatum)
    idDatum = createDatumplane(modelName, partName[0], datumPlane = YZPLANE, datumDepth = 0.0)
    createPartion(modelName, partName[0], idDatum)
    ExtrudPartition(modelName, partName[0])
    #create pile geometry
    createPart_Pile(modelName, partName, pileDia, pileLength)
    idDatum = createDatumplane(modelName, partName[1], datumPlane = YZPLANE, datumDepth = 0.0)
    createPartion(modelName, partName[1], idDatum)
    #define and assign material properties 
    createPSet_cellWhole(modelName, partName[0], 'setSoil')
    defineProperties_soil(modelName, materialName, sectionName, Es0, EsL, nu_s, pileLength)
    assignProperties_soil(modelName, partName[0], sectionName, 'setSoil')
    createPSet_cellWhole(modelName, partName[1], 'setPile')
    defineProperties_pile(modelName, materialName, sectionName, Ep, nu_p)
    assignProperties_pile(modelName, partName[1], sectionName, 'setPile')
    #Assembly, Setps and Output
    createAssembly(modelName, partName, instanceName)
    createStep(modelName)
    createOutputReq(modelName, fieldoutputReqName)
    # Interaction
    # Tie constraint b/n pile and soil
    surfRegion_s1 = createISurface_region(modelName, instanceName[0], 'mSurf_1', (pileDia/2.0)*cos(3.0/4.0*pi), -(pileDia/2.0)*sin(3.0/4.0*pi), pileLength/2.0)
    surfRegion_s2 = createISurface_region(modelName, instanceName[0], 'mSurf_2', (pileDia/2.0)*cos(1.0/4.0*pi), -(pileDia/2.0)*sin(1.0/4.0*pi), pileLength/2.0)
    surfRegion_s3 = createISurface_region(modelName, instanceName[0], 'mSurf_3', ((pileDia/2.0)*cos(3.0/4.0*pi))/2.0, -((pileDia/2.0)*sin(3.0/4.0*pi))/2.0, pileLength)
    surfRegion_s4 = createISurface_region(modelName, instanceName[0], 'mSurf_4', ((pileDia/2.0)*cos(1.0/4.0*pi))/2.0, -((pileDia/2.0)*sin(1.0/4.0*pi))/2.0, pileLength)
    surfRegion_p1 = createISurface_region(modelName, instanceName[1], 'sSurf_1', ((pileDia/2.0)*cos(3.0/4.0*pi)), -((pileDia/2.0)*sin(3.0/4.0*pi)), pileLength/2.0)
    surfRegion_p2 = createISurface_region(modelName, instanceName[1], 'sSurf_2', ((pileDia/2.0)*cos(1.0/4.0*pi)), -((pileDia/2.0)*sin(1.0/4.0*pi)), pileLength/2.0)
    surfRegion_p3 = createISurface_region(modelName, instanceName[1], 'sSurf_3', ((pileDia/2.0)*cos(3.0/4.0*pi))/2.0, -((pileDia/2.0)*sin(3.0/4.0*pi))/2.0, pileLength)
    surfRegion_p4 = createISurface_region(modelName, instanceName[1], 'sSurf_4', ((pileDia/2.0)*cos(1.0/4.0*pi))/2.0, -((pileDia/2.0)*sin(1.0/4.0*pi))/2.0, pileLength)
    createConstraint_tie(modelName, surfRegion_p1, surfRegion_s1, 'Tie-1')
    createConstraint_tie(modelName, surfRegion_p2, surfRegion_s2, 'Tie-2')
    # createConstraint_tie(modelName, surfRegion_s3, surfRegion_p3, 'Tie-3')
    # createConstraint_tie(modelName, surfRegion_s4, surfRegion_p4, 'Tie-4') #REMOVED as BASE TIED TO THE SOIL EXIHIBTED SLIGHTLY HIGHER STIFNESS COMPARED TO BVP FOR SHORT RIGID PILES
    # definition of sets
    # TODO: find a better method to define and used sets or to locate various regions of the model
    #sets on the base
    createPSet_face(modelName, partName[0], 'baseFace_1', ((domainDia/2.0)*cos(3.0/4.0*pi))/2.0, -((domainDia/2.0)*sin(3.0/4.0*pi))/2.0, domainDepth)
    createPSet_face(modelName, partName[0], 'baseFace_2', ((domainDia/2.0)*cos(1.0/4.0*pi))/2.0, -((domainDia/2.0)*sin(1.0/4.0*pi))/2.0, domainDepth)
    createPSet_face(modelName, partName[0], 'mbaseFace_1', ((pileDia/2.0)*cos(3.0/4.0*pi))/2.0, -((pileDia/2.0)*sin(3.0/4.0*pi))/2.0, domainDepth)
    createPSet_face(modelName, partName[0], 'mbaseFace_2', ((pileDia/2.0)*cos(1.0/4.0*pi))/2.0, -((pileDia/2.0)*sin(1.0/4.0*pi))/2.0, domainDepth)
    createPSet_edge(modelName, partName[0], 'baseEdge_1', (domainDia/2.0)*cos(3.0/4.0*pi), -(domainDia/2.0)*sin(3.0/4.0*pi), domainDepth)
    createPSet_edge(modelName, partName[0], 'baseEdge_2', (domainDia/2.0)*cos(1.0/4.0*pi), -(domainDia/2.0)*sin(1.0/4.0*pi), domainDepth)
    createPSet_edge(modelName, partName[0], 'mbaseEdge_1', (pileDia/2.0)*cos(3.0/4.0*pi), -(pileDia/2.0)*sin(3.0/4.0*pi), domainDepth)
    createPSet_edge(modelName, partName[0], 'mbaseEdge_2', (pileDia/2.0)*cos(1.0/4.0*pi), -(pileDia/2.0)*sin(1.0/4.0*pi), domainDepth)
    createPSet_edge(modelName, partName[0], 'baseEdge_3', 0.0, -domainDia/4.0, domainDepth)
    createPSet_edge(modelName, partName[0], 'mbaseEdge_3', 0.0, -pileDia/4.0, domainDepth)
    createPSet_edge(modelName, partName[0], 'baseEdge_4', -domainDia/4.0, 0.0, domainDepth)
    createPSet_edge(modelName, partName[0], 'baseEdge_5', domainDia/4.0, 0.0, domainDepth)
    createPSet_edge(modelName, partName[0], 'mbaseEdge_4', -pileDia/4.0, 0.0, domainDepth)
    createPSet_edge(modelName, partName[0], 'mbaseEdge_5', pileDia/4.0, 0.0, domainDepth)
    #sets on the top
    createPSet_face(modelName, partName[0], 'topFace_1', ((domainDia/2.0)*cos(3.0/4.0*pi))/2.0, -((domainDia/2.0)*sin(3.0/4.0*pi))/2.0, 0.0)
    createPSet_face(modelName, partName[0], 'topFace_2', ((domainDia/2.0)*cos(1.0/4.0*pi))/2.0, -((domainDia/2.0)*sin(1.0/4.0*pi))/2.0, 0.0)
    createPSet_edge(modelName, partName[0], 'topEdge_1', (domainDia/2.0)*cos(3.0/4.0*pi), -(domainDia/2.0)*sin(3.0/4.0*pi), 0.0)
    createPSet_edge(modelName, partName[0], 'topEdge_2', (domainDia/2.0)*cos(1.0/4.0*pi), -(domainDia/2.0)*sin(1.0/4.0*pi), 0.0)
    createPSet_edge(modelName, partName[0], 'mtopEdge_1', (pileDia/2.0)*cos(3.0/4.0*pi), -(pileDia/2.0)*sin(3.0/4.0*pi), 0.0)
    createPSet_edge(modelName, partName[0], 'mtopEdge_2', (pileDia/2.0)*cos(1.0/4.0*pi), -(pileDia/2.0)*sin(1.0/4.0*pi), 0.0)
    createPSet_edge(modelName, partName[0], 'topEdge_3', 0.0, -domainDia/4.0, 0.0)
    createPSet_edge(modelName, partName[0], 'topEdge_4', -domainDia/4.0, 0.0, 0.0)
    createPSet_edge(modelName, partName[0], 'topEdge_5', domainDia/4.0, 0.0, 0.0)
    #sets on front
    createPSet_face(modelName, partName[0], 'frontFace_1', -(domainDia/2.0)/2.0, 0.0, pileLength/2.0)
    createPSet_face(modelName, partName[0], 'frontFace_2', (domainDia/2.0)/2.0, 0.0, pileLength/2.0)
    createPSet_face(modelName, partName[0], 'frontFace_3', -(domainDia/2.0)/2.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    createPSet_face(modelName, partName[0], 'frontFace_4', (domainDia/2.0)/2.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    #mfrontFace_1 = mSurf_1
    #mfrontFace_2 = mSurf_2
    createPSet_face(modelName, partName[0], 'mfrontFace_1', (pileDia/2.0)*cos(3.0/4.0*pi), -(pileDia/2.0)*sin(3.0/4.0*pi), pileLength/2.0)
    createPSet_face(modelName, partName[0], 'mfrontFace_2', (pileDia/2.0)*cos(1.0/4.0*pi), -(pileDia/2.0)*sin(1.0/4.0*pi), pileLength/2.0)
    createPSet_face(modelName, partName[0], 'mfrontFace_3', -(pileDia/2.0)/2.0, -(pileDia/2.0)/2.0, pileLength)
    createPSet_face(modelName, partName[0], 'mfrontFace_4', (pileDia/2.0)/2.0, -(pileDia/2.0)/2.0, pileLength)
    createPSet_face(modelName, partName[0], 'mfrontFace_5', -(pileDia/2.0)/2.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    createPSet_face(modelName, partName[0], 'mfrontFace_6', (pileDia/2.0)/2.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    createPSet_edge(modelName, partName[0], 'frontEdge_1', domainDia/2.0, 0.0, pileLength/2.0)
    createPSet_edge(modelName, partName[0], 'frontEdge_2', pileDia/2.0, 0.0, pileLength/2.0)
    createPSet_edge(modelName, partName[0], 'frontEdge_3', -pileDia/2.0, 0.0, pileLength/2.0)
    createPSet_edge(modelName, partName[0], 'frontEdge_4', -domainDia/2.0, 0.0, pileLength/2.0)
    createPSet_edge(modelName, partName[0], 'frontEdge_5', domainDia/2.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    createPSet_edge(modelName, partName[0], 'frontEdge_6', pileDia/2.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    createPSet_edge(modelName, partName[0], 'frontEdge_7', -pileDia/2.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    createPSet_edge(modelName, partName[0], 'frontEdge_8', -domainDia/2.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    createPSet_edge(modelName, partName[0], 'frontEdge_9', -domainDia/4.0, 0.0, pileLength)
    createPSet_edge(modelName, partName[0], 'frontEdge_10', domainDia/4.0, 0.0, pileLength)
    createPSet_edge(modelName, partName[0], 'mfrontEdge_1', 0.0, -pileDia/2.0, pileLength/2.0)
    createPSet_edge(modelName, partName[0], 'mfrontEdge_4', (pileDia/2.0)*cos(3.0/4.0*pi), -(pileDia/2.0)*sin(3.0/4.0*pi), pileLength)
    createPSet_edge(modelName, partName[0], 'mfrontEdge_5', (pileDia/2.0)*cos(1.0/4.0*pi), -(pileDia/2.0)*sin(1.0/4.0*pi), pileLength)
    createPSet_edge(modelName, partName[0], 'mfrontEdge_6', -pileDia/4.0, 0.0, pileLength)
    createPSet_edge(modelName, partName[0], 'mfrontEdge_7', pileDia/4.0, 0.0, pileLength)
    createPSet_edge(modelName, partName[0], 'mfrontEdge_3', 0.0, -pileDia/4.0, pileLength)
    createPSet_edge(modelName, partName[0], 'mfrontEdge_2', 0.0, 0.0, pileLength + (domainDepth-pileLength)/2.0/2.0)
    #sets on back
    createPSet_face(modelName, partName[0], 'backFace_1', (domainDia/2.0)*cos(3.0/4.0*pi), -(domainDia/2.0)*sin(3.0/4.0*pi), pileLength/2.0)
    createPSet_face(modelName, partName[0], 'backFace_2', (domainDia/2.0)*cos(1.0/4.0*pi), -(domainDia/2.0)*sin(1.0/4.0*pi), pileLength/2.0)
    createPSet_face(modelName, partName[0], 'backFace_3', (domainDia/2.0)*cos(3.0/4.0*pi), -(domainDia/2.0)*sin(3.0/4.0*pi), pileLength + (domainDepth-pileLength)/2.0)
    createPSet_face(modelName, partName[0], 'backFace_4', (domainDia/2.0)*cos(1.0/4.0*pi), -(domainDia/2.0)*sin(1.0/4.0*pi), pileLength + (domainDepth-pileLength)/2.0)
    createPSet_edge(modelName, partName[0], 'backEdge_1', 0.0, -domainDia/2.0, pileLength/2.0)
    createPSet_edge(modelName, partName[0], 'backEdge_2', 0.0, -domainDia/2.0, pileLength + (domainDepth-pileLength)/2.0)
    createPSet_edge(modelName, partName[0], 'backEdge_3', (domainDia/2.0)*cos(1.0/4.0*pi), -(domainDia/2.0)*sin(3.0/4.0*pi), pileLength)
    createPSet_edge(modelName, partName[0], 'backEdge_4', (domainDia/2.0)*cos(3.0/4.0*pi), -(domainDia/2.0)*sin(3.0/4.0*pi), pileLength)
    #sets on pile
    createPSet_face(modelName, partName[1], 'pfrontFace_1', -pileDia/4.0, 0.0, pileLength/2.0)
    createPSet_face(modelName, partName[1], 'pfrontFace_2', pileDia/4.0, 0.0, pileLength/2.0)
    createPSet_edge(modelName, partName[1], 'ptopEdge_1', -pileDia/4.0, 0.0, 0.0)
    createPSet_edge(modelName, partName[1], 'ptopEdge_2', pileDia/4.0, 0.0, 0.0)
    createPSet_edge(modelName, partName[1], 'pbaseEdge_1', -pileDia/4.0, 0.0, pileLength)
    createPSet_edge(modelName, partName[1], 'pbaseEdge_2', pileDia/4.0, 0.0, pileLength)
    createPSet_edge(modelName, partName[1], 'pfrontEdge_1', -pileDia/2.0, 0.0, pileLength/2.0)
    createPSet_edge(modelName, partName[1], 'pfrontEdge_2', pileDia/2.0, 0.0, pileLength/2.0)
    createPSet_edge(modelName, partName[1], 'pfrontEdge_3', 0.0, 0.0, pileLength/2.0)
    #Fix the base of the domain
    baseFxsets = ['baseFace_1', 'baseFace_2', 'baseEdge_4', 'baseEdge_5', 'mbaseEdge_4', 'mbaseEdge_5', 'mbaseFace_1', 'mbaseFace_2', 'baseEdge_1', 'baseEdge_2', 'mbaseEdge_1', 'mbaseEdge_2', 'baseEdge_3', 'mbaseEdge_3', 'baseEdge_4', 'baseEdge_5', 'mbaseEdge_4', 'mbaseEdge_5']
    for m in range(len(baseFxsets)):
        assignBC_BaseFx(modelName, instanceName[0], baseFxsets[m], baseFxsets[m])
    #Fix the lateral sides
    latFxsets = ['backFace_1', 'backFace_2', 'backFace_3', 'backFace_4', 'topEdge_1', 'topEdge_2', 'backEdge_1', 'backEdge_2', 'backEdge_3', 'backEdge_4']
    for n in range(len(latFxsets)):
        assignBC_lateralFx(modelName, instanceName[0], latFxsets[n], latFxsets[n])
    #Apply Y-symmetry 
    YsymmSets_1 = ['frontFace_1', 'frontFace_2', 'frontFace_3', 'frontFace_4', 'frontEdge_1', 'frontEdge_2', 'frontEdge_3', 'frontEdge_4', 'frontEdge_5', 'frontEdge_6', 'frontEdge_7', 'frontEdge_8', 'frontEdge_9', 'frontEdge_10', 'mfrontEdge_2', 'mfrontEdge_6', 'mfrontEdge_7', 'topEdge_4', 'topEdge_5', 'mfrontFace_5', 'mfrontFace_6']
    for o in range(len(YsymmSets_1)):
        assignBC_Ysymm(modelName, instanceName[0], YsymmSets_1[o], YsymmSets_1[o])
    YsymmSets_2 = ['pfrontFace_1', 'pfrontFace_2', 'ptopEdge_1', 'ptopEdge_2', 'pbaseEdge_1', 'pbaseEdge_2', 'pfrontEdge_1', 'pfrontEdge_2', 'pfrontEdge_3']
    for r in range(len(YsymmSets_2)):
        assignBC_Ysymm(modelName, instanceName[1], YsymmSets_2[r], YsymmSets_2[r])
    #apply the pile head
    createISurface_region(modelName, instanceName[1], 'pileHead_1', ((pileDia/2.0)*cos(3.0/4.0*pi))/2.0, -((pileDia/2.0)*sin(3.0/4.0*pi))/2.0, 0.0)
    createISurface_region(modelName, instanceName[1], 'pileHead_2', ((pileDia/2.0)*cos(1.0/4.0*pi))/2.0, -((pileDia/2.0)*sin(3.0/4.0*pi))/2.0, 0.0)
    # # Create Constraint to for moment application
    # refPoint1, refPoint_Loc1 = CreateRefPoint(modelName, 'pileHead_RP', 0.0, 0.0, -0.5)
    # createConstraint_Coupling(modelName, 'pileHead_RP', 'pileHead_1', 'Coupling-1')
    # createConstraint_Coupling(modelName, 'pileHead_RP', 'pileHead_2', 'Coupling-2')
    # #
    # assignPileHeadMoment(modelName, 'pileHead_RP', 'PileHeadMoment_1', pileHeadMoment)
    #
    assignPileHeadLatLoad(modelName, instanceName[1], 'pileHead_1', 'pileHeadLoad_1', pileHeadTraction)
    assignPileHeadLatLoad(modelName, instanceName[1], 'pileHead_2', 'pileHeadLoad_2', pileHeadTraction)
    #assign dummy temp field
    assingDTempField(modelName, instanceName[0], 'setSoil')
    #assign mesh controls
    assignMeshControlWEDGE(modelName, instanceName[0], -pileDia/4.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    assignMeshControlWEDGE(modelName, instanceName[0], pileDia/4.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    assignMeshControlWEDGE(modelName, instanceName[1], -pileDia/4.0, 0.0, pileLength/2.0)
    assignMeshControlWEDGE(modelName, instanceName[1], pileDia/4.0, 0.0, pileLength/2.0)
    #seed the cylinderical edges of the soil
    createSeed_EdgeNum(modelName, instanceName[0], seedNumber, (domainDia/2.0)*cos(3.0/4.0*pi), -(domainDia/2.0)*sin(3.0/4.0*pi), 0.0)
    createSeed_EdgeNum(modelName, instanceName[0], seedNumber, (domainDia/2.0)*cos(1.0/4.0*pi), -(domainDia/2.0)*sin(1.0/4.0*pi), 0.0)
    createSeed_EdgeNum(modelName, instanceName[0], seedNumber, (pileDia/2.0)*cos(3.0/4.0*pi), -(pileDia/2.0)*sin(3.0/4.0*pi), 0.0)
    createSeed_EdgeNum(modelName, instanceName[0], seedNumber, (pileDia/2.0)*cos(1.0/4.0*pi), -(pileDia/2.0)*sin(1.0/4.0*pi), 0.0)
    createSeed_EdgeNum(modelName, instanceName[0], seedNumber, (domainDia/2.0)*cos(1.0/4.0*pi), -(domainDia/2.0)*sin(3.0/4.0*pi), pileLength)
    createSeed_EdgeNum(modelName, instanceName[0], seedNumber, (domainDia/2.0)*cos(3.0/4.0*pi), -(domainDia/2.0)*sin(3.0/4.0*pi), pileLength)
    createSeed_EdgeNum(modelName, instanceName[0], seedNumber, (pileDia/2.0)*cos(3.0/4.0*pi), -(pileDia/2.0)*sin(3.0/4.0*pi), pileLength)
    createSeed_EdgeNum(modelName, instanceName[0], seedNumber, (pileDia/2.0)*cos(1.0/4.0*pi), -(pileDia/2.0)*sin(1.0/4.0*pi), pileLength)
    createSeed_EdgeNum(modelName, instanceName[0], seedNumber, (domainDia/2.0)*cos(1.0/4.0*pi), -(domainDia/2.0)*sin(3.0/4.0*pi), domainDepth)
    createSeed_EdgeNum(modelName, instanceName[0], seedNumber, (domainDia/2.0)*cos(3.0/4.0*pi), -(domainDia/2.0)*sin(3.0/4.0*pi), domainDepth)
    createSeed_EdgeNum(modelName, instanceName[0], seedNumber, (pileDia/2.0)*cos(3.0/4.0*pi), -(pileDia/2.0)*sin(3.0/4.0*pi), domainDepth)
    createSeed_EdgeNum(modelName, instanceName[0], seedNumber, (pileDia/2.0)*cos(1.0/4.0*pi), -(pileDia/2.0)*sin(1.0/4.0*pi), domainDepth)
    #seed the radial edges of soil
    createSeed_EdgeBias_2(modelName, instanceName[0], seedSizeMin_R, seedSizeMax_R, 0.0, -domainDia/4.0, 0.0)
    createSeed_EdgeBias_2(modelName, instanceName[0], seedSizeMin_R, seedSizeMax_R, -domainDia/4.0, 0.0, 0.0)
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_R, seedSizeMax_R, domainDia/4.0, 0.0, 0.0)
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_R, seedSizeMax_R, 0.0, -domainDia/4.0, pileLength)
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_R, seedSizeMax_R, -domainDia/4.0, 0.0, pileLength)
    createSeed_EdgeBias_2(modelName, instanceName[0], seedSizeMin_R, seedSizeMax_R, domainDia/4.0, 0.0, pileLength)
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_R, seedSizeMax_R, 0.0, -domainDia/4.0, domainDepth)
    createSeed_EdgeBias_2(modelName, instanceName[0], seedSizeMin_R, seedSizeMax_R, -domainDia/4.0, 0.0, domainDepth)
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_R, seedSizeMax_R, domainDia/4.0, 0.0, domainDepth)
    #sed the vertical edges below the pile
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_Z, seedSizeMax_Z, domainDia/2.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_Z, seedSizeMax_Z, pileDia/2.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    createSeed_EdgeBias_2(modelName, instanceName[0], seedSizeMin_Z, seedSizeMax_Z, -pileDia/2.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_Z, seedSizeMax_Z, -domainDia/2.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_Z, seedSizeMax_Z, 0.0, 0.0, pileLength + (domainDepth-pileLength)/2.0)
    createSeed_EdgeBias_2(modelName, instanceName[0], seedSizeMin_Z, seedSizeMax_Z, 0.0, -domainDia/2.0, pileLength + (domainDepth-pileLength)/2.0)
    #seed the radial edges of semi-circle regions under the pile
    createSeed_EdgeNum(modelName, instanceName[0], 1, -pileDia/4.0, 0.0, pileLength)
    createSeed_EdgeNum(modelName, instanceName[0], 1, pileDia/4.0, 0.0, pileLength)
    createSeed_EdgeNum(modelName, instanceName[0], 1, 0.0, -pileDia/4.0, pileLength)
    createSeed_EdgeNum(modelName, instanceName[0], 1, -pileDia/4.0, 0.0, domainDepth)
    createSeed_EdgeNum(modelName, instanceName[0], 1, pileDia/4.0, 0.0, domainDepth)
    createSeed_EdgeNum(modelName, instanceName[0], 1, 0.0, -pileDia/4.0, domainDepth)
    #seed the perimeter of the top and bottom of pile
    createSeed_EdgeNum(modelName, instanceName[1], seedNumber, (pileDia/2.0)*cos(1.0/4.0*pi), -(pileDia/2.0)*sin(3.0/4.0*pi), 0.0)
    createSeed_EdgeNum(modelName, instanceName[1], seedNumber, (pileDia/2.0)*cos(3.0/4.0*pi), -(pileDia/2.0)*sin(3.0/4.0*pi), 0.0)
    createSeed_EdgeNum(modelName, instanceName[1], seedNumber, (pileDia/2.0)*cos(1.0/4.0*pi), -(pileDia/2.0)*sin(3.0/4.0*pi), pileLength)
    createSeed_EdgeNum(modelName, instanceName[1], seedNumber, (pileDia/2.0)*cos(3.0/4.0*pi), -(pileDia/2.0)*sin(3.0/4.0*pi), pileLength)
    #seed the radial edges of the pile top and bottom
    createSeed_EdgeNum(modelName, instanceName[1], 1, -pileDia/4.0, 0.0, pileLength)
    createSeed_EdgeNum(modelName, instanceName[1], 1, pileDia/4.0, 0.0, pileLength)
    createSeed_EdgeNum(modelName, instanceName[1], 1, 0.0, -pileDia/4.0, pileLength)
    createSeed_EdgeNum(modelName, instanceName[1], 1, -pileDia/4.0, 0.0, 0.0)
    createSeed_EdgeNum(modelName, instanceName[1], 1, pileDia/4.0, 0.0, 0.0)
    createSeed_EdgeNum(modelName, instanceName[1], 1, 0.0, -pileDia/4.0, 0.0)
    # apply global seed
    #createSeed_global(modelName, instanceName, seedSizeGlob)
    # seed along the length of the pile and soil
    createSeed_EdgeBias_2(modelName, instanceName[1], seedSizeMin_Z_1, seedSizeMax_Z_1, 0.0, 0.0, pileLength/2.0)
    createSeed_EdgeBias_2(modelName, instanceName[1], seedSizeMin_Z_1, seedSizeMax_Z_1, -pileDia/2.0, 0.0, pileLength/2.0)
    createSeed_EdgeBias_2(modelName, instanceName[1], seedSizeMin_Z_1, seedSizeMax_Z_1, pileDia/2.0, 0.0, pileLength/2.0)
    createSeed_EdgeBias_2(modelName, instanceName[1], seedSizeMin_Z_1, seedSizeMax_Z_1, 0.0, -pileDia/2.0, pileLength/2.0)
    #
    createSeed_EdgeBias_2(modelName, instanceName[0], seedSizeMin_Z_1, seedSizeMax_Z_1, -pileDia/2.0, 0.0, pileLength/2.0)
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_Z_1, seedSizeMax_Z_1, pileDia/2.0, 0.0, pileLength/2.0)
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_Z_1, seedSizeMax_Z_1, 0.0, -pileDia/2.0, pileLength/2.0)
    #
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_Z_1, seedSizeMax_Z_1, -domainDia/2.0, 0.0, pileLength/2.0)
    createSeed_EdgeBias_1(modelName, instanceName[0], seedSizeMin_Z_1, seedSizeMax_Z_1, domainDia/2.0, 0.0, pileLength/2.0)
    createSeed_EdgeBias_2(modelName, instanceName[0], seedSizeMin_Z_1, seedSizeMax_Z_1, 0.0, -domainDia/2.0, pileLength/2.0)
    #
    # create mesh
    createMesh(modelName, instanceName)
    #assign elements to soil
    #hex elements
    assignElemType_Hex(modelName, instanceName[0], ((domainDia/2.0)*cos(3.0/4.0*pi))/2.0, -((domainDia/2.0)*sin(3.0/4.0*pi))/2.0, pileLength/2.0)
    assignElemType_Hex(modelName, instanceName[0], ((domainDia/2.0)*cos(1.0/4.0*pi))/2.0, -((domainDia/2.0)*sin(1.0/4.0*pi))/2.0, pileLength/2.0)
    assignElemType_Hex(modelName, instanceName[0], ((domainDia/2.0)*cos(3.0/4.0*pi))/2.0, -((domainDia/2.0)*sin(3.0/4.0*pi))/2.0, pileLength + (domainDepth-pileLength)/2.0)
    assignElemType_Hex(modelName, instanceName[0], ((domainDia/2.0)*cos(1.0/4.0*pi))/2.0, -((domainDia/2.0)*sin(1.0/4.0*pi))/2.0, pileLength + (domainDepth-pileLength)/2.0)
    #wedge elements
    assignElemType_Wedge(modelName, instanceName[0], ((pileDia/2.0)*cos(3.0/4.0*pi))/2.0, -((pileDia/2.0)*sin(3.0/4.0*pi))/2.0, pileLength + (domainDepth-pileLength)/2.0)
    assignElemType_Wedge(modelName, instanceName[0], ((pileDia/2.0)*cos(1.0/4.0*pi))/2.0, -((pileDia/2.0)*sin(1.0/4.0*pi))/2.0, pileLength + (domainDepth-pileLength)/2.0)
    #assign elements to Pile
    #only wedge elements
    assignElemType_Wedge(modelName, instanceName[1], ((pileDia/2.0)*cos(3.0/4.0*pi))/2.0, -((pileDia/2.0)*sin(3.0/4.0*pi))/2.0, pileLength/2.0)
    assignElemType_Wedge(modelName, instanceName[1], ((pileDia/2.0)*cos(1.0/4.0*pi))/2.0, -((pileDia/2.0)*sin(1.0/4.0*pi))/2.0, pileLength/2.0)
    #create job
    createJob(jobName, modelName)
    #submit job
    submitJob(jobName)

