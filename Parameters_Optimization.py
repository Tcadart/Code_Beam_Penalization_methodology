# Version python 3.9.7
# ----------------------------------------------------------------------------
# Created By  : Thomas Cadart
# Created Date: 17/04/2023
# version ='1.0'
# ---------------------------------------------------------------------------
""" Code for studying lattice structures based on density and simulation method.
This code has highlighted the issues with finite element simulation of Timoshenko beams
and developed a penalty methodology to enable the use of beam elements for the study of
lattice structures across a wide range of densities and geometries.
This code works with Abaqus to carry out finite element simulations.
To compare the lattice structures, homogenization was necessary.
To simplify its implementation, the Micromechanics plugin was used.
Link: https://www.linkedin.com/pulse/micromechanics-plugin-abaqus-ross-mclendon/
"""
# ---------------------------------------------------------------------------
import sys
import os
current_directory = os.getcwd()
sys.path.append(current_directory)
from Lattice_description import *
from Material import *

from abaqus import *
from abaqusConstants import *
import regionToolset
import xyPlot
from operator import add

sys.path.insert(8, r"D:/travail_Abaqus/MicroMechanics_v1.18/MicroMechanics")
from microMechanics.mmpBackend import Interface
from microMechanics.mmpBackend.mmpInterface.mmpRVEConstants import *
from microMechanics.mmpBackend.mmpKernel.mmpLibrary import *
from microMechanics.mmpBackend import mmpKernel as Kernel
from odbAccess import openOdb

import numpy as np
from math import *
from numpy import linalg as LA
from scipy.optimize import minimize_scalar
import re
# *******************************************************************************************************************
# *******************************************************************************************************************

#                                           Functions

# *******************************************************************************************************************
# *******************************************************************************************************************

def Create_Struct(Geometry,Cell_dim,Cell_number,name_model,Frequency_analysis,Maillage_size=0.5):
    Kernel.Library.Lattice.generateGeneralLattice(name_model,(Cell_dim,Cell_dim,Cell_dim),Geometry,(Cell_number, Cell_number, Cell_number),Maillage_size)

def Load_And_Job(name_Job,name_model,Frequency_analysis=0):
    #Create a job
    mdb.Job(name=name_Job, model=name_model, description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
    multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
    #Frequency 
    if Frequency_analysis == 1:
        #Loading
        output = (False, False, False)
        defPeriodicBC(name_model,name_Job,output)
        #Frequency step
        mdb.models[name_model].FrequencyStep(name='DummyStep', previous='Initial', 
        maintainAttributes=True, numEigen=20)
        mdb.models[name_model].boundaryConditions['PinnedNode'].suppress()
        # Boundary condition
        a = mdb.models[name_model].rootAssembly
        region = a.sets['RP-Normal']
        mdb.models[name_model].DisplacementBC(name='Periodic_Normal', createStepName='DummyStep', 
            region=region, u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)
        region = a.sets['RP-Shear']
        mdb.models[name_model].DisplacementBC(name='Periodic_Shear', createStepName='DummyStep', 
            region=region, u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)
        #Job
        mdb.jobs[name_Job].submit(consistencyChecking=OFF)
        mdb.jobs[name_Job].waitForCompletion()
    if Frequency_analysis == 0:
        #Assign loads and submit job
        output = (True, False, True)
        defPeriodicBC(name_model,name_Job,output)
        #Submit and wait
        job1 = mdb.Job(name=name_Job,model=name_model)
        job1.submit()
        job1.waitForCompletion()

def defPeriodicBC(name_model,name_Job,output):
    Interface.Loading.MechanicalModelMaker(constraintType='PERIODIC', 
    drivenField='STRAIN', modelName=name_model, jobName=name_Job, 
    doNotSubmit=True, homogenizeProperties=output)

def Post_Process_homogenization(name_Job,name_model):
    current_directory = os.getcwd()
    Interface.PostProcess.MechanicalPostProcessWorkflow(model=name_model, 
    ODBName=current_directory+"/"+str(name_Job)+".odb", doHomogenization=True,
    materialType=('Engineering Constants', ), fieldAveraging=False, 
    getStrainConcentration=False, averageVolume=False, 
    averageVolumeBySection=False, selectedFields=(), 
    getStatisticalDistribution=False)

def Get_matrix(name_model):
    mat = list(mdb.models[name_model].materials['Homog_EngConst_t=0pt0'].elastic.table)
    return mat

def Homogenization_3D_Solid(Lattice_geom,Cell_number,name_model,name_Job,Radius,Frequency_analysis,Maillage_size=0.5):
    # Create Model
    Create_Struct(Lattice_geom,1.0,Cell_number,name_model,Frequency_analysis,Maillage_size)
    #Create Load and Job and Submit for analysis
    Load_And_Job(name_Job,name_model,Frequency_analysis)
    if Frequency_analysis == 1:
        return 0
    else:
        # Create Homogenization matrix
        Post_Process_homogenization(name_Job,name_model)
        #Get homogenization matrix
        mat = Get_matrix(name_model)
        vol = Get_Volume(name_model)
        result = np.array([Radius,vol,mat[0][0],mat[0][1],mat[0][2],mat[0][3],mat[0][4],mat[0][5],mat[0][6],mat[0][7],mat[0][8]])
        return result

def Homogenization_3D_Beam(Lattice_geom,name_model,name_Job,Radius,name_Part,VectorOrientation,CorrectionExteriorBeam,Frequency_analysis,mesh_size,material_type):
    CreateModel(name_model)
    CreatePart(name_model,name_Part)
    Create_geometry_beam(Lattice_geom,name_model,name_Part)
    mesh_beam(name_model,name_Part,mesh_size)
    Section_beam(name_model,name_Part,Radius,VectorOrientation,CorrectionExteriorBeam,Frequency_analysis,material_type)
    Assembly_beam(name_model,name_Part,'LatticeCell-1')
    Load_And_Job(name_Job,name_model,Frequency_analysis)
    if Frequency_analysis == 0:
        Post_Process_homogenization(name_Job,name_model)
        mat = Get_matrix(name_model)
        vol = Get_Volume(name_model)
        result = np.array([Radius,vol,mat[0][0],mat[0][1],mat[0][2],mat[0][3],mat[0][4],mat[0][5],mat[0][6],mat[0][7],mat[0][8]])
        return result
    else:
        return 0

def Homogenization_3D_BeamMod(Lattice_geom,name_model,name_Job,Radius,name_Part,PourcentMod,RadiusFactor,VectorOrientation,CorrectionExteriorBeam,Frequency_analysis,material_type,mesh_size=0.05,MassSolid=0):
    CreateModel(name_model)
    CreatePart(name_model,name_Part)
    Create_geometry_beamMod(Lattice_geom,name_model,PourcentMod,name_Part)
    mesh_beam(name_model,name_Part,mesh_size)
    Section_beam_Mod(name_model,name_Part,Radius,RadiusFactor,Lattice_geom,VectorOrientation,CorrectionExteriorBeam,Frequency_analysis,material_type)
    Assembly_beam(name_model,name_Part,'LatticeCell-1')
    if Frequency_analysis == 1 and MassSolid != 0: # Case density mod frequency
        ChangeDensityBeamMod(name_model,MassSolid,material_type)
    Load_And_Job(name_Job,name_model,Frequency_analysis)
    if Frequency_analysis == 0:
        Post_Process_homogenization(name_Job,name_model)
        mat = Get_matrix(name_model)
        vol = Get_Volume(name_model)
        result = np.array([Radius,vol,mat[0][0],mat[0][1],mat[0][2],mat[0][3],mat[0][4],mat[0][5],mat[0][6],mat[0][7],mat[0][8]],np.mean(PourcentMod))
        return result
    else:
        return 0

def Create_geometry_beam(Lattice_geom,name_model,name_Part):
    p = mdb.models[name_model].parts[name_Part]
    Wire = Lattice_geom[0]
    List_point = []
    List_point.append(Wire[0:3])
    List_beam = np.zeros((len(Lattice_geom),2))
    for i in range (len(Lattice_geom)):
        counter_1 = 0
        counter_2 = 0
        Wire = Lattice_geom[i]
        for j in range(len(List_point)):
            if (Wire[0:3] == List_point[j][0:3]):
                counter_1 = counter_1+1
                List_beam[i][0] = j
            if (Wire[3:6] == List_point[j][0:3]):
                counter_2 = counter_2+1
                List_beam[i][1] = j
            if counter_1==1 and counter_2==1:
                break
            if (j == len(List_point)-1):
                if counter_1 == 0:
                    List_point.append(Wire[0:3])
                    List_beam[i][0] = j+1
                if counter_2 == 0:
                    List_point.append(Wire[3:6])
                    if counter_1 == 1:
                        List_beam[i][1] = j+1
                    else:
                        List_beam[i][1] = j+2
                break 
    #Create each point
    for i in range(len(List_point)):
        p.DatumPointByCoordinate(coords=(List_point[i][0], List_point[i][1], List_point[i][2]))
    #Create each beam
    d2 = p.datums
    for i in range(len(List_beam)):
        a = 2 + int(List_beam[i][0])
        b = 2 + int(List_beam[i][1])
        if (CorrectionExteriorBeam == 2) and (((List_point[int(List_beam[i][0])][0]==List_point[int(List_beam[i][1])][0]) and List_point[int(List_beam[i][0])][0]>=1) or ((List_point[int(List_beam[i][0])][1]==List_point[int(List_beam[i][1])][1]) and List_point[int(List_beam[i][0])][1]>=1) or ((List_point[int(List_beam[i][0])][2]==List_point[int(List_beam[i][1])][2]) and List_point[int(List_beam[i][0])][2]>=1)):
            a = 0
        else:
            p.WirePolyLine(points=((d2[a], d2[b]), ), mergeType=IMPRINT, meshable=ON)

def mesh_beam(name_model,name_Part,mesh_size=0.05):
    p = mdb.models[name_model].parts[name_Part]
    p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()


def Section_beam(name_model,name_Part,Radius,VectorOrientation,CorrectionExteriorBeam,Frequency_analysis,material_type):
    #Beam Section Orientation
    BeamSectionOrientation(name_model,name_Part,VectorOrientation)
    # Beam Profile creation
    BeamProfileDefinition(name_model,Radius,CorrectionExteriorBeam,0,RadiusFactor)
    #Create Material
    name_material = create_material(material_type,name_model,Frequency_analysis)
    # Beam Section creation
    mdb.models[name_model].BeamSection(name='CircBeams', integration=DURING_ANALYSIS, 
    poissonRatio=0.0, profile='Circ', material=name_material,
    temperatureVar=LINEAR, consistentMassMatrix=False)
    if CorrectionExteriorBeam == 1:
        mdb.models[name_model].BeamSection(name='CircBeamsExt', integration=DURING_ANALYSIS, 
        poissonRatio=0.0, profile='CircExt', material=name_material,
        temperatureVar=LINEAR, consistentMassMatrix=False)
    # Beam Section assignment
    if CorrectionExteriorBeam != 1:
        p = mdb.models[name_model].parts[name_Part]
        e = p.edges
        edges = e.getSequenceFromMask(mask=('[#ffffffff:7 #3fff ]',),)
        region = p.Set(edges=edges, name='AllBeams')
        p = mdb.models[name_model].parts[name_Part]
        p.SectionAssignment(region=region, sectionName='CircBeams', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', 
            thicknessAssignment=FROM_SECTION)
    elif CorrectionExteriorBeam == 1:
        # Beam Section assignment with correction
        p = mdb.models[name_model].parts[name_Part]
        e=p.edges
        edges_ext = []
        for i in range(len(Lattice_geom)):
            #Check if beam is exterior
            if (((Lattice_geom[i][0]==Lattice_geom[i][3]) and Lattice_geom[i][0]==1) or ((Lattice_geom[i][1]==Lattice_geom[i][4]) and Lattice_geom[i][1]==1) or ((Lattice_geom[i][2]==Lattice_geom[i][5]) and Lattice_geom[i][2]==1)):
                edges_ext.append(e.findAt(((Lattice_geom[i][0]+(Lattice_geom[i][3]-Lattice_geom[i][0])/2,Lattice_geom[i][1]+(Lattice_geom[i][4]-Lattice_geom[i][1])/2,Lattice_geom[i][2]+(Lattice_geom[i][5]-Lattice_geom[i][2])/2),),))
        region_ext = p.Set(edges=edges_ext, name='BeamExt')
        p.SectionAssignment(region=region_ext, sectionName='CircBeamsExt', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', 
            thicknessAssignment=FROM_SECTION)
        # For Interior beams
        p = mdb.models[name_model].parts[name_Part]
        region_int = p.SetByBoolean(name = 'BeamInt',operation=DIFFERENCE,sets=(p.sets['AllBeams'],p.sets['BeamExt'],))
        p.SectionAssignment(region=region_int, sectionName='CircBeams', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', 
            thicknessAssignment=FROM_SECTION)

def Section_beam_Mod(name_model,name_Part,Radius,RadiusFactor,Lattice_geom,VectorOrientation,CorrectionExteriorBeam,Frequency_analysis,material_type):
    #Beam Section Orientation
    BeamSectionOrientation(name_model,name_Part,VectorOrientation)
    # Beam Profile creation
    BeamProfileDefinition(name_model,Radius,CorrectionExteriorBeam,1,RadiusFactor)
    #Create Material
    name_material = create_material(material_type,name_model,Frequency_analysis)
    # Beam Section creation
    mdb.models[name_model].BeamSection(name='CircBeams', integration=DURING_ANALYSIS,
    poissonRatio=0.0, profile='Circ', material=name_material,
    temperatureVar=LINEAR, consistentMassMatrix=False)
    mdb.models[name_model].BeamSection(name='CircBeamsMod', integration=DURING_ANALYSIS,
    poissonRatio=0.0, profile='Circ_Mod', material=name_material,
    temperatureVar=LINEAR, consistentMassMatrix=False)
    if CorrectionExteriorBeam == 1:
        mdb.models[name_model].BeamSection(name='CircBeamsExt', integration=DURING_ANALYSIS,
        poissonRatio=0.0, profile='CircExt', material=name_material,
        temperatureVar=LINEAR, consistentMassMatrix=False)
    # Beam Section assignment
    if CorrectionExteriorBeam !=1:
        # For normal beams
        p = mdb.models[name_model].parts[name_Part]
        e=p.edges
        edges_mid = []
        for i in range(len(Lattice_geom)):
            edges_mid.append(e.findAt(((Lattice_geom[i][0]+(Lattice_geom[i][3]-Lattice_geom[i][0])/2,Lattice_geom[i][1]+(Lattice_geom[i][4]-Lattice_geom[i][1])/2,Lattice_geom[i][2]+(Lattice_geom[i][5]-Lattice_geom[i][2])/2),),))
        region_mid = p.Set(edges=edges_mid, name='BeamMid')
        p.SectionAssignment(region=region_mid, sectionName='CircBeams', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
        # For Modified beams
        p = mdb.models[name_model].parts[name_Part]
        region_ext = p.SetByBoolean(name = 'BeamMod',operation=DIFFERENCE,sets=(p.sets['AllBeams'],p.sets['BeamMid'],))
        p.SectionAssignment(region=region_ext, sectionName='CircBeamsMod', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
    elif CorrectionExteriorBeam == 1:
        p = mdb.models[name_model].parts[name_Part]
        e=p.edges
        edges_mid = []
        edges_ext = []
        for i in range(len(Lattice_geom)):
            #Check if beam is exterior
            if (((Lattice_geom[i][0]==Lattice_geom[i][3]) and Lattice_geom[i][0]==1) or ((Lattice_geom[i][1]==Lattice_geom[i][4]) and Lattice_geom[i][1]==1) or ((Lattice_geom[i][2]==Lattice_geom[i][5]) and Lattice_geom[i][2]==1)):
                edges_ext.append(e.findAt(((Lattice_geom[i][0]+(Lattice_geom[i][3]-Lattice_geom[i][0])/2,Lattice_geom[i][1]+(Lattice_geom[i][4]-Lattice_geom[i][1])/2,Lattice_geom[i][2]+(Lattice_geom[i][5]-Lattice_geom[i][2])/2),),))
                edges_ext.append(e.findAt(((Lattice_geom[i][0]+(Lattice_geom[i][3]-Lattice_geom[i][0])/100,Lattice_geom[i][1]+(Lattice_geom[i][4]-Lattice_geom[i][1])/100,Lattice_geom[i][2]+(Lattice_geom[i][5]-Lattice_geom[i][2])/100),),))
                edges_ext.append(e.findAt(((Lattice_geom[i][3]-(Lattice_geom[i][3]-Lattice_geom[i][0])/100,Lattice_geom[i][4]-(Lattice_geom[i][4]-Lattice_geom[i][1])/100,Lattice_geom[i][5]-(Lattice_geom[i][5]-Lattice_geom[i][2])/100),),))
            else:
                edges_mid.append(e.findAt(((Lattice_geom[i][0]+(Lattice_geom[i][3]-Lattice_geom[i][0])/2,Lattice_geom[i][1]+(Lattice_geom[i][4]-Lattice_geom[i][1])/2,Lattice_geom[i][2]+(Lattice_geom[i][5]-Lattice_geom[i][2])/2),),))
        # For Exterior Beams
        region_ext = p.Set(edges=edges_ext, name='BeamExt')
        p.SectionAssignment(region=region_ext, sectionName='CircBeamsExt', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
        # For Middle Beams
        region_mid = p.Set(edges=edges_mid, name='BeamMid')
        p.SectionAssignment(region=region_mid, sectionName='CircBeams', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
        # For Modified Beams
        p = mdb.models[name_model].parts[name_Part]
        region_ext = p.SetByBoolean(name = 'BeamMod',operation=DIFFERENCE,sets=(p.sets['AllBeams'],p.sets['BeamMid'],p.sets['BeamExt'],))
        p.SectionAssignment(region=region_ext, sectionName='CircBeamsMod', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

def BeamSectionOrientation(name_model,name_Part,VectorOrientation):
    p = mdb.models[name_model].parts[name_Part]
    e = p.edges
    edges = e.getSequenceFromMask(mask=('[#ffffffff:7 #3fff ]',),)
    region = p.Set(edges=edges, name='AllBeams')
    region=regionToolset.Region(edges=edges)
    p = mdb.models[name_model].parts[name_Part]
    p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(VectorOrientation[0], VectorOrientation[1], VectorOrientation[2]))

def BeamProfileDefinition(name_model,Radius,CorrectionExteriorBeam,BeamMod,RadiusFactor=0):
    mdb.models[name_model].CircularProfile(name='Circ', r=Radius)
    if BeamMod == 1:
        mdb.models[name_model].CircularProfile(name='Circ_Mod', r=Radius*RadiusFactor)
    if CorrectionExteriorBeam == 1:
        mdb.models[name_model].CircularProfile(name='CircExt', r=Radius/50) #Create model for exterior beams

def Assembly_beam(name_model,name_Part,name_Assembly):
    a1 = mdb.models[name_model].rootAssembly
    a1.DatumCsysByDefault(CARTESIAN)
    p = mdb.models[name_model].parts[name_Part]
    a1.Instance(name=name_Assembly, part=p, dependent=ON)


def CreateModel(name_model):
    #Create Model
    mdb.Model(name=name_model, modelType=STANDARD_EXPLICIT)

def CreatePart(name_model,name_Part):
    #Create Part
    s1 = mdb.models[name_model].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.Line(point1=(0.0, 0.0), point2=(0.0, 5.0))
    s1.VerticalConstraint(entity=g[2], addUndoState=False)
    p = mdb.models[name_model].Part(name=name_Part, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models[name_model].parts[name_Part]
    p.BaseWire(sketch=s1)
    s1.unsetPrimaryObject()
    p = mdb.models[name_model].parts[name_Part]
    del mdb.models[name_model].sketches['__profile__']
    # Delete Part
    del p.features['Wire-1']

def Get_Volume(name_model):
    vol=mdb.models[name_model].rootAssembly.getMassProperties()
    return vol['volume']

def Get_Volume_region(name_model, set_name):
    a = mdb.models[name_model].rootAssembly
    partSet = a.sets[set_name]
    edge_array = partSet.edges
    region = regionToolset.Region(edges=edge_array)
    massprop = a.getMassProperties(region)
    return massprop['volume']

def Get_Mass_region(name_model, set_name):
    a = mdb.models[name_model].rootAssembly
    partSet = a.sets[set_name]
    edge_array = partSet.edges
    region = regionToolset.Region(edges=edge_array)
    massprop = a.getMassProperties(region)
    return massprop['mass']

def ChangeDensityBeamMod(name_model,MassSolid,material_type):
    MassBeamMid = Get_Mass_region(name_model,'LatticeCell-1.BeamMid')
    VolBeamMod = Get_Volume_region(name_model,'LatticeCell-1.BeamMod')
    DensiteMod = (float(MassSolid)-float(MassBeamMid))/float(VolBeamMod)
    name_material = create_material(material_type,name_model,20,DensiteMod)
    mdb.models[name_model].sections['CircBeamsMod'].setValues(poissonRatio=0.0, material=name_material)


def Create_geometry_beamMod(Lattice_geom,name_model,PourcentMod,name_Part):
    p = mdb.models[name_model].parts[name_Part]
    Wire = Lattice_geom[0]
    List_point = []
    List_point_Mod = []
    Position_point_Mod = np.zeros((len(Lattice_geom),6))
    List_point.append(Wire[0:3])
    List_beam = np.zeros((len(Lattice_geom),2))
    for i in range (len(Lattice_geom)):
        counter_1 = 0
        counter_2 = 0
        Wire = Lattice_geom[i]
        for j in range(len(List_point)):
            #Detect if point already added
            if (Wire[0:3] == List_point[j][0:3]):
                counter_1 = counter_1+1
                List_beam[i][0] = j
            if (Wire[3:6] == List_point[j][0:3]):
                counter_2 = counter_2+1
                List_beam[i][1] = j
            if counter_1==1 and counter_2==1:
                break
            #Add new point
            if (j == len(List_point)-1):
                if counter_1 == 0:
                    List_point.append(Wire[0:3])
                    List_beam[i][0] = j+1
                if counter_2 == 0:
                    List_point.append(Wire[3:6])
                    if counter_1 == 1:
                        List_beam[i][1] = j+1
                    else:
                        List_beam[i][1] = j+2
                break
            #Calculate intermediate point
        Position_point_Mod[i] = IntermediatePoint(Wire,PourcentMod,i)
    #Create each point
    nPoint = 0
    for i in range(len(List_point)):
        p.DatumPointByCoordinate(coords=(List_point[i][0], List_point[i][1], List_point[i][2]))
        nPoint = nPoint + 1
    for i in range(len(Position_point_Mod)):
        nPoint = nPoint + 1
        List_point_Mod.append([nPoint,nPoint+1])
        p.DatumPointByCoordinate(coords=(Position_point_Mod[i][0], Position_point_Mod[i][1], Position_point_Mod[i][2]))
        nPoint = nPoint + 1
        p.DatumPointByCoordinate(coords=(Position_point_Mod[i][3], Position_point_Mod[i][4], Position_point_Mod[i][5]))
    #Create each beam
    d2 = p.datums
    for i in range(len(List_beam)):
        if (CorrectionExteriorBeam == 2) and (((List_point[int(List_beam[i][0])][0]==List_point[int(List_beam[i][1])][0]) and List_point[int(List_beam[i][0])][0]>=1) or ((List_point[int(List_beam[i][0])][1]==List_point[int(List_beam[i][1])][1]) and List_point[int(List_beam[i][0])][1]>=1) or ((List_point[int(List_beam[i][0])][2]==List_point[int(List_beam[i][1])][2]) and List_point[int(List_beam[i][0])][2]>=1)):
            a = 0
        else:
            valCorr = 2
            a = valCorr + int(List_beam[i][0])
            b = valCorr + int(List_point_Mod[i][0]) - 1
            p.WirePolyLine(points=((d2[a], d2[b]), ), mergeType=IMPRINT, meshable=ON)
            a = valCorr + int(List_point_Mod[i][0]) - 1
            b = valCorr + int(List_point_Mod[i][1]) - 1
            p.WirePolyLine(points=((d2[a], d2[b]), ), mergeType=IMPRINT, meshable=ON)
            a = valCorr + int(List_point_Mod[i][1]) - 1
            b = valCorr + int(List_beam[i][1])
            p.WirePolyLine(points=((d2[a], d2[b]), ), mergeType=IMPRINT, meshable=ON)

def IntermediatePoint(Wire,PourcentMod,i):
    DR =[Wire[3]-Wire[0],Wire[4]-Wire[1],Wire[5]-Wire[2]]
    Position_point_Mod = np.zeros(6)
    for j in range(3):
        Position_point_Mod[j] = Wire[j]+DR[j]*PourcentMod[i][0]/100
        Position_point_Mod[j+3] = Wire[j]+DR[j]*(1.0-PourcentMod[i][1]/100.0)
    return Position_point_Mod


def GetPourcentMod(Radius_beam,Lattice_geom,TypePourcentMod):
    LengthModtemp = np.zeros((len(Lattice_geom),2))
    LengthMod = GetmatlengthMod(Lattice_geom)
    for i in range(len(Lattice_geom)):
        for j in range(2):
            if LengthMod[i][j] <=0.01:
                LengthModtemp[i][j] = 0.01
            else:
                LengthBeam = sqrt((Lattice_geom[i][3] - Lattice_geom[i][0]) ** 2 + (Lattice_geom[i][4] - Lattice_geom[i][1]) ** 2 + (Lattice_geom[i][5] - Lattice_geom[i][2]) ** 2)
                if TypePourcentMod == 'Constant':
                    LengthModtemp[i][j] = 10
                if TypePourcentMod == 'Variable':
                    LengthModtemp[i][j] = round((Radius_beam/LengthBeam)*100,2)
                if TypePourcentMod == 'Variable2':
                    LengthModtemp[i][j] = round((sqrt(3)*(Radius_beam/LengthBeam))*100,2)
                if TypePourcentMod == 'Variablecut':
                    LengthModtemp[i][j] = LengthMod[i][j]*100
    return LengthModtemp


def GetRelativeError(HomoSolid,HomoBeamMod):
    data_process_s = np.array([1/HomoSolid[2],1/HomoSolid[3],1/HomoSolid[4],1/HomoSolid[8],1/HomoSolid[9],1/HomoSolid[10],-HomoSolid[5]/HomoSolid[2],-HomoSolid[5]/HomoSolid[3],-HomoSolid[6]/HomoSolid[2],-HomoSolid[6]/HomoSolid[4],-HomoSolid[7]/HomoSolid[3],-HomoSolid[7]/HomoSolid[4]])
    data_process_b = np.array([1/HomoBeamMod[2],1/HomoBeamMod[3],1/HomoBeamMod[4],1/HomoBeamMod[8],1/HomoBeamMod[9],1/HomoBeamMod[10],-HomoBeamMod[5]/HomoBeamMod[2],-HomoBeamMod[5]/HomoBeamMod[3],-HomoBeamMod[6]/HomoBeamMod[2],-HomoBeamMod[6]/HomoBeamMod[4],-HomoBeamMod[7]/HomoBeamMod[3],-HomoBeamMod[7]/HomoBeamMod[4]])
    mat_rigid = np.array(
        [[data_process_s[0],data_process_s[6],data_process_s[8], 0, 0, 0],
        [data_process_s[7],data_process_s[1],data_process_s[10], 0, 0, 0],
        [data_process_s[9],data_process_s[11],data_process_s[2], 0, 0, 0],
        [ 0, 0, 0, data_process_s[3], 0, 0],
        [ 0, 0, 0, 0, data_process_s[4], 0],
        [ 0, 0, 0, 0, 0, data_process_s[5]]]
    )
    mat_rigid_inv_s = LA.inv(mat_rigid)
    normfrob_s = LA.norm(mat_rigid_inv_s,'fro')
    mat_rigid = np.array(
        [[data_process_b[0],data_process_b[6],data_process_b[8], 0, 0, 0],
        [data_process_b[7],data_process_b[1],data_process_b[10], 0, 0, 0],
        [data_process_b[9],data_process_b[11],data_process_b[2], 0, 0, 0],
        [ 0, 0, 0, data_process_b[3], 0, 0],
        [ 0, 0, 0, 0, data_process_b[4], 0],
        [ 0, 0, 0, 0, 0, data_process_b[5]]]
    )
    mat_rigid_inv_b = LA.inv(mat_rigid)
    #Calculation norm dif
    mat_rigid_dif = mat_rigid_inv_s - mat_rigid_inv_b
    normfrob_dif = LA.norm(mat_rigid_dif,'fro')
    #Calculation relative error
    Erreur_relative = normfrob_dif/normfrob_s
    return Erreur_relative


def FuncOptimize(RadiusFactor,HomoSolid,Lattice_geom,name_model,name_Job,Radius_beam,name_Part,PourcentMod,VectorOrientation,CorrectionExteriorBeam,material_type):
    if isinstance(RadiusFactor, np.ndarray):
        RadiusFactor = RadiusFactor[0]
    HomoBeamMod = Homogenization_3D_BeamMod(Lattice_geom,name_model,name_Job,Radius_beam,name_Part,PourcentMod,RadiusFactor,VectorOrientation,CorrectionExteriorBeam,0,material_type)
    Relative_Error = np.array([GetRelativeError(HomoSolid,HomoBeamMod)])
    return Relative_Error


def Getangle(liaison,Lattice_geom):
    angle = []
    angle_deg = []
    for j in range(len(liaison)):
        u = [Lattice_geom[liaison[0]][3]-Lattice_geom[liaison[0]][0],Lattice_geom[liaison[0]][4]-Lattice_geom[liaison[0]][1],Lattice_geom[liaison[0]][5]-Lattice_geom[liaison[0]][2]]
        v = [Lattice_geom[liaison[j]][3]-Lattice_geom[liaison[j]][0],Lattice_geom[liaison[j]][4]-Lattice_geom[liaison[j]][1],Lattice_geom[liaison[j]][5]-Lattice_geom[liaison[j]][2]]
        if np.dot(u, v) < (np.linalg.norm(u) * np.linalg.norm(v)):
            angle_rad = acos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
            angle_deg.append(degrees(angle_rad))
        else:
            angle_deg.append(0)
    angle_deg = np.array(angle_deg)
    if np.all(angle_deg == 0):
        angle.append(0)
    else:
        non_zero_angle = [x for x in angle_deg if x >= 0.01]
        angle.append(min(non_zero_angle))
    angle = np.array(angle)
    return angle

def GetmatlengthMod(Lattice_geom):
    #Step 1 Find node in lattice
    point = []
    for i in range(len(Lattice_geom)):
        if list(Lattice_geom[i])[0:3] not in point:
            point.append(list(Lattice_geom[i])[0:3])
        if list(Lattice_geom[i])[3:6] not in point:
            point.append(list(Lattice_geom[i])[3:6])
    # Find Min and max point
    minpoint = 0.5
    maxpoint = 0.5
    cornerpoint = []
    exteriorpointX = []
    exteriorpointY = []
    exteriorpointZ = []
    for i in range(len(Lattice_geom)):
        for j in range(6):
            if Lattice_geom[i][j]<minpoint:
                minpoint = Lattice_geom[i][j]
            if Lattice_geom[i][j]>maxpoint:
                maxpoint = Lattice_geom[i][j]
    # Find corner point and exterior
    for i in range(len(point)):
        if (point[i][0] == minpoint or point[i][0] == maxpoint) and (point[i][1] == minpoint or point[i][1] == maxpoint) and (point[i][2] == minpoint or point[i][2] == maxpoint):
            cornerpoint.append(point[i])
        elif point[i][0] == minpoint or point[i][0] == maxpoint:
            exteriorpointX.append(point[i])
        elif point[i][1] == minpoint or point[i][1] == maxpoint:
            exteriorpointY.append(point[i])
        elif point[i][2] == minpoint or point[i][2] == maxpoint:
            exteriorpointZ.append(point[i])
    liaison = []
    liaisoncorner = []
    liaisonexteriorX = []
    liaisonexteriorY = []
    liaisonexteriorZ = []
    for i in range(len(point)):
        liaison.append([])
        for j in range(len(Lattice_geom)):
            if point[i] == list(Lattice_geom[j])[0:3] or point[i] == list(Lattice_geom[j])[3:6]:
                if point[i] in cornerpoint:
                    liaisoncorner.append(j)
                elif point[i] in exteriorpointX:
                    liaisonexteriorX.append(j)
                elif point[i] in exteriorpointY:
                    liaisonexteriorY.append(j)
                elif point[i] in exteriorpointZ:
                    liaisonexteriorZ.append(j)
                else:
                    liaison[i].append(j)
    cornerpoint = np.array(cornerpoint)
    angle = []
    angle_deg = []
    for i in range(len(liaison)):
        angle.append(0)
        if len(liaison[i])>1:
            for j in range(len(liaison[i])):
                    u = [Lattice_geom[i][3]-Lattice_geom[i][0],Lattice_geom[i][4]-Lattice_geom[i][1],Lattice_geom[i][5]-Lattice_geom[i][2]]
                    v = [Lattice_geom[liaison[i][j]][3]-Lattice_geom[liaison[i][j]][0],Lattice_geom[liaison[i][j]][4]-Lattice_geom[liaison[i][j]][1],Lattice_geom[liaison[i][j]][5]-Lattice_geom[liaison[i][j]][2]]
                    if np.dot(u, v) < (np.linalg.norm(u) * np.linalg.norm(v)):
                        angle_rad = acos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
                        angle_deg.append(degrees(angle_rad))
                    else:
                        angle_deg.append(0)
            non_zero_angle = [x for x in angle_deg if x != 0]
            angle[i] = min(non_zero_angle)
    if len(liaisoncorner)>1:
        anglecorner = Getangle(liaisoncorner,Lattice_geom)
    if len(liaisonexteriorX)>1:
        angleexteriorX = Getangle(liaisonexteriorX,Lattice_geom)
    if len(liaisonexteriorY)>1:
        angleexteriorY = Getangle(liaisonexteriorY,Lattice_geom)
    if len(liaisonexteriorZ)>1:
        angleexteriorZ = Getangle(liaisonexteriorZ,Lattice_geom)
    pointangle = np.zeros(len(point))
    for i in range(len(point)):
        if any(np.array_equal(point[i], elem) for elem in cornerpoint):
            pointangle[i] = anglecorner
        elif any(np.array_equal(point[i], elem) for elem in exteriorpointX):
            pointangle[i] = angleexteriorX
        elif any(np.array_equal(point[i], elem) for elem in exteriorpointY):
            pointangle[i] = angleexteriorY
        elif any(np.array_equal(point[i], elem) for elem in exteriorpointZ):
            pointangle[i] = angleexteriorZ
        elif len(liaison[i])!=0:
            pointangle[i] = angle[i]
    LengthMod = GetLengthModWithAngle(pointangle,Lattice_geom,point)
    return LengthMod

def GetLengthModWithAngle(pointangle,Lattice_geom,point):
    LengthModtemp = np.zeros(len(pointangle))
    for i in range(len(pointangle)):
        if pointangle[i]!=0 and pointangle[i]!=180:
            LengthModtemp[i] = Lattice_geom[0][-1] / tan((pointangle[i] / 2) * pi / 180)
        else:
            LengthModtemp[i]=0.00001
    LengthMod = np.zeros((len(Lattice_geom),2))
    for i in range(len(Lattice_geom)):
        for j in range(len(point)):
            if tuple(point[j]) ==  Lattice_geom[i][0:3]:
                LengthMod[i][0] = LengthModtemp[j]
            if tuple(point[j]) ==  Lattice_geom[i][3:6]:
                LengthMod[i][1] = LengthModtemp[j]
    return LengthMod


def Visualization_Frequency(name_Job):
    current_directory = os.getcwd()
    session.viewports['Viewport: 1'].setValues(
        displayedObject=session.odbs[current_directory+'/'+name_Job+'.odb'])
    o3 = session.openOdb(
        name=current_directory+'/'+name_Job+'.odb')
    session.viewports['Viewport: 1'].setValues(displayedObject=o3)
    session.viewports['Viewport: 1'].makeCurrent()


def getFrequencyResult(name_Job):
    current_directory = os.getcwd()
    odb = openOdb(path=current_directory+'/'+name_Job+'.odb')
    xy0 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Eigenfrequency: EIGFREQ for Whole Model', 
    suppressQuery=True, __linkedVpName__='Viewport: 1')
    freq = session.Curve(xyData=xy0)
    xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Eigenvalue: EIGVAL for Whole Model', 
    suppressQuery=True, __linkedVpName__='Viewport: 1')
    val = session.Curve(xyData=xy1)
    freqData = freq.data
    valData = val.data
    result = [(x[0], x[1], y[1]) for x, y in zip(valData, freqData)]
    odb.close()
    return result 

def delete_all_models():
    model_names = list(mdb.models.keys())
    for name in model_names:
        if name != 'Model-1':  # 'Model-1' basic model
            del mdb.models[name]



def LimitConditionPlasticity(name_model,name_Job, limit_condition):
    defaut = ['Default' if limit_condition[i] != '' else '' for i in range(6)]
    Interface.Loading.MechanicalModelMaker(constraintType='PERIODIC', 
    drivenField='STRAIN', modelName=name_model, 
    mechanicalHistoryType=LOADUSERDEFINED, fieldHistory1=(str(limit_condition[0]), str(limit_condition[1]), str(limit_condition[2]), str(limit_condition[3]), 
    str(limit_condition[4]), str(limit_condition[5]), defaut[0], defaut[1], defaut[2], defaut[3], defaut[4], defaut[5]), homogenizationIntervals=0, 
    jobName=name_Job, doNotSubmit=True, homogenizeProperties=(False, False, 
    False), totalHistoryTime=1)

def delUnnecessary(name_model):
    del mdb.models[name_model].fieldOutputRequests['F-Output-1']
    del mdb.models[name_model].historyOutputRequests['H-Output-1']


def SimuPlasticity(Lattice_geom,Cell_dim, Cell_number,name_model,name_Job,Radius_beam,Frequency_analysis,Maillage_size,NLgeom,limit_condition,TimeInc,materialType):
    Create_Struct(Lattice_geom, 1.0, Cell_number, name_model, Frequency_analysis, Maillage_size*10)
    # Create Job
    mdb.Job(name=name_Job, model=name_model, description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
        multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
    # Material
    name_material = create_material(materialType, name_model, Frequency_analysis)
    mdb.models[name_model].sections['Ti-6Al-4V'].setValues(material=name_material,
                                                          thickness=None)
    # Periodic Boundary Condition
    LimitConditionPlasticity(name_model,name_Job, limit_condition)
    name_step = 'LoadHist-1'
    if NLgeom == 1:
        mdb.models[name_model].steps[name_step].setValues(nlgeom=ON)
    timeIncrement(name_model, name_step,TimeInc)
    delUnnecessary(name_model)
    # Launch Job
    Submit_Job(name_Job,name_model)
    # get result
    vol = Get_Volume(name_model)
    dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime = getPlasticityData(name_Job,name_step)
    return dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime,vol

def Submit_Job(name_Job,name_model):
    #Submit and wait
    job1 = mdb.Job(name=name_Job,model=name_model)
    job1.submit()
    job1.waitForCompletion()

def getPlasticityData(name_Job, name_step):
    # Open the ODB file and extract data
    odb = openOdb(name_Job + '.odb')
    historyRegionKey = list(odb.steps[name_step].historyRegions.keys())[0]
    historyRegion = odb.steps[name_step].historyRegions[historyRegionKey]
    data_keys = ['RF1', 'RF2', 'RF3', 'U1', 'U2', 'U3']
    results = {key: [item[1] for item in historyRegion.historyOutputs[key].data] for key in data_keys}
    dataTime = [item[0] for item in historyRegion.historyOutputs['U3'].data]
    odb.close()
    # Extract CPU Time from the .msg file
    with open(name_Job + '.msg', 'r') as file:
        for last_line in file:
            pass
    CPUTime = float(re.findall(r'\d+\.\d+|\d+', last_line)[0])
    return_tuple = tuple(results[key] for key in data_keys) + (dataTime, CPUTime)
    return return_tuple

def PlasticityBeam(name_model,name_Part,Lattice_geom,mesh_size,Radius,VectorOrientation,CorrectionExteriorBeam,Frequency_analysis,name_Job,Type,NLgeom,limit_condition,TimeInc):
    CreateModel(name_model)
    CreatePart(name_model,name_Part)
    # Geometry lattice
    if Type == 'Beam':
        Create_geometry_beam(Lattice_geom,name_model,name_Part)
    elif Type == 'BeamMod' or Type == 'BeamModVariable':
        Create_geometry_beamMod(Lattice_geom,name_model,PourcentMod,name_Part)
    mesh_beam(name_model,name_Part,mesh_size)
    # Section Beam
    if Type == 'Beam':
        Section_beam(name_model,name_Part,Radius,VectorOrientation,CorrectionExteriorBeam,1,material_type)
    elif Type == 'BeamMod' or Type == 'BeamModVariable':
        Section_beam_Mod(name_model,name_Part,Radius,RadiusFactor,Lattice_geom,VectorOrientation,CorrectionExteriorBeam,Frequency_analysis,material_type)
    # Assembly lattice
    name_Assembly = 'LatticeCell-1'
    Assembly_beam(name_model,name_Part,name_Assembly)
    # Create Job
    mdb.Job(name=name_Job, model=name_model, description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
        multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
    # Periodic Boundary Condition
    LimitConditionPlasticity(name_model,name_Job,limit_condition)
    name_step = 'LoadHist-1'
    if NLgeom == 1:
        mdb.models[name_model].steps[name_step].setValues(nlgeom=ON)
    timeIncrement(name_model, name_step,TimeInc)
    changeBeamOutput(name_model)
    delUnnecessary(name_model)
    # Launch Job
    Submit_Job(name_Job,name_model)
    dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime = getPlasticityData(name_Job,name_step)
    return dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime

def timeIncrement(name_model,name_step,TimeInc = 0.05):
    mdb.models[name_model].steps[name_step].setValues(
        timeIncrementationMethod=FIXED, initialInc=TimeInc, noStop=OFF)

def changeBeamOutput(name_model):
    mdb.models[name_model].fieldOutputRequests['LoadHistory_FieldOutput'].setValues(
        variables=('S', 'E', 'PE', 'NE', 'LE', 'U', 'RF', 'CF', 'SF', 'TRSHR',
                   'TRNOR', 'TEMP', 'IVOL'))

def WriteResult(Radius,dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime):
    f.writelines([str(Radius) + "\n"])
    f.writelines([str(CPUTime) + "\n"])
    data_arrays = [dataTime, dataRF1, dataRF2, dataRF3, dataU1, dataU2, dataU3]
    for i, data in enumerate(data_arrays, start=1):
        f.writelines([np.array2string(np.array(data), separator=',') + "\n"])


def determine_file_name(lattice_type, analysis_type, is_plasticity, is_nlgeom, is_frequency_analysis):

    file_name = Type_lattice(lattice_type) + "_" + analysis_type

    if analysis_type == 'FrequencyMod':
        file_name += 'Dens_Freq'
    else:
        if is_plasticity:
            file_name += "Plasticity"
            if is_nlgeom:
                file_name += "NL"
        if is_frequency_analysis:
            file_name += "Freq"

    file_name += ".txt"

    return file_name

def handle_visualization(name_model, n, visual):
    if visual == 1:
        return name_model +'_'+ str(n)
    else:
        del mdb.models[name_model]
        del mdb.jobs[name_job]

def write_results_to_file(f, lattice_type, name_job, radius_beam, result, is_plasticity, frequency_analysis):
    if frequency_analysis:
        frequency_result = getFrequencyResult(name_job)
        vol = Get_Volume(name_model)
        f.writelines([str(vol)+"\n"])
        for item in frequency_result:
            f.write(','.join(map(str, item)) + '\n')
    elif is_plasticity:
        dataRF1, dataRF2, dataRF3, dataU1, dataU2, dataU3, dataTime, CPUTime, vol = result
        WriteResult(vol, dataRF1, dataRF2, dataRF3, dataU1, dataU2, dataU3, dataTime, CPUTime)
    else:
        if radius_beam is not None:
            result_str = ','.join(map(str, [radius_beam] + list(result)))
            f.writelines([result_str + "\n"])
        else:
            f.writelines([np.array2string(result, separator=',')+"\n"])


def perform_analysis(type_analysis, lattice_geom, name_model, name_job, radius_beam, name_part, pourcent_mod,
                     radius_factor, vector_orientation, correction_exterior_beam, mesh_size, frequency_analysis, nlgeom,
                     limit_condition, time_inc, cell_number, cell_dim, is_plasticity,material_type):
    if is_plasticity:
        if type_analysis == 'Solid':
            return SimuPlasticity(lattice_geom, cell_dim, cell_number, name_model, name_job, radius_beam,
                                  frequency_analysis, mesh_size, nlgeom, limit_condition, time_inc,material_type)
        # Add similar conditionals for other types like 'Beam', 'BeamMod', etc. if needed
        return PlasticityBeam(name_model, name_part, lattice_geom, mesh_size, radius_beam, vector_orientation,
                              correction_exterior_beam, 1, name_job, type_analysis, nlgeom, limit_condition, time_inc)

    elif type_analysis == 'Solid':
        return Homogenization_3D_Solid(lattice_geom, cell_number, name_model, name_job, radius_beam, frequency_analysis,
                                       mesh_size)

    elif type_analysis == 'Beam':
        return Homogenization_3D_Beam(lattice_geom, name_model, name_job, radius_beam, name_part, vector_orientation,
                                      correction_exterior_beam, frequency_analysis, mesh_size,material_type)

    elif type_analysis in ['BeamMod', 'BeamModVariable']:
        return Homogenization_3D_BeamMod(lattice_geom, name_model, name_job, radius_beam, name_part, pourcent_mod,
                                         radius_factor, vector_orientation, correction_exterior_beam,material_type,
                                         frequency_analysis, mesh_size)



#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Variables

#*******************************************************************************************************************
#*******************************************************************************************************************

Radius_beam = np.array([0.015,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2])
# Radius_beam = np.array([0.05]) # For test
Cell_number = 1
name_model = 'Lattice'
name_Job = 'Job_1'
name_Part = 'LatticeCell'
Lattice = 0
# 0 => BCC 
# 1 => Octet 
# 2 => OctetExt 
# 3 => OctetInt 
# 4 => BCCZ 
# 5 => Cubic 
# 6 => OctahedronZ 
# 7 => OctahedronZcross 
# 8 => Kelvin 
# 9 => Cubic formulation 2 (centered)
# 10 => Cubic V3
# 11 => Cubic V4
# 12 => New lattice
# 13 => Diamond

Type = 'Solid' #Solid or Beam or BeamMod or BeamModVariable or Optimization_bounded or FrequencyMod
TypePourcentMod = 'Variablecut' #Type of calculation of modified length
#Constant : 10 %
#Variable : Radius
#Variable2 : sqrt(3)*Radius
#Variablecut : Calculation with angle at node

RadiusFactor = 1.5 #Penalization coefficient
#Visualization lattice
Visual = 1 # 0 : Inactive / 1 : Active
deletevar = 0
TypeCorrection = 2
# Correction = 1 / Reduction exterior beams
# Correction = 2 / Delete exterior beams
Frequency_analysis = 0
# 1 -> active
Plasticity = 1
# 1 = On / 0 = Off
# Non-linear behavior
NLgeom = 1
# 1 = On / 0 = Off
# Mesh size
mesh_size = 0.05
# Displacement direction
# [E11,E22,E33,E23,E13,E12]
limit_condition = ['','','-0.3','','','']
TimeInc = 0.025
material_type = 2
#0 : Ti-6Al-4V
#1 : VeroClear

delete_all_models()
CorrectionExteriorBeam = GetCorrectionExteriorBeam(Lattice,TypeCorrection)
VectorOrientation = getVectorOrientation(Lattice)
file_name = determine_file_name(Lattice, Type, Plasticity, NLgeom, Frequency_analysis)
f = open(file_name, "w")

if (Type in ['Solid', 'Beam', 'BeamMod', 'BeamModVariable']):
    for n in range(len(Radius_beam)):
        Lattice_geom = Lattice_geometry_corrected(Lattice, Radius_beam[n])
        PourcentMod = GetPourcentMod(Radius_beam[n], Lattice_geom, TypePourcentMod) if Type in ['BeamMod','BeamModVariable'] else None
        result = perform_analysis(Type, Lattice_geom, name_model, name_Job, Radius_beam[n], name_Part, PourcentMod,
                                  RadiusFactor, VectorOrientation, CorrectionExteriorBeam, mesh_size,
                                  Frequency_analysis, NLgeom, limit_condition, TimeInc, Cell_number, 1.0,
                                  Plasticity,material_type)
        write_results_to_file(f, Lattice, name_Job, Radius_beam[n], result, Plasticity, Frequency_analysis)
        handle_visualization(name_model, n, Visual)

###################### Frequency Modifier of beam density ######################
if (Type == 'FrequencyMod'):
    Type == 'BeamModVariable'
    for n in range(len(Radius_beam)):
        # Get Volume of solid model
        Lattice_geom = Lattice_geometry_corrected(Lattice, Radius_beam[n])
        name_model_Solid = name_model+'Solid'
        Create_Struct(Lattice_geom,1.0,Cell_number,name_model_Solid,Frequency_analysis,Maillage_size)
        volSolid = Get_Volume(name_model_Solid)
        MassSolid = volSolid*4.429e-09
        # Create BeamMod Model
        Lattice_geom = Lattice_geometry_corrected(Lattice, Radius_beam[n])
        PourcentMod = GetPourcentMod(Radius_beam[n],Lattice_geom,TypePourcentMod)
        Homogenization_3D_BeamMod(Lattice_geom,name_model,name_Job,Radius_beam[n],name_Part,PourcentMod,RadiusFactor,VectorOrientation,CorrectionExteriorBeam,Frequency_analysis,material_type,mesh_size,MassSolid)
        result = getFrequencyResult(name_Job)
        write_results_to_file(f, Lattice, name_Job, Radius_beam[n], result, Plasticity,Frequency_analysis)
        handle_visualization(name_model, n, Visual)

###################### Parameters Optimization Scipy ###################### 
if Type == 'Optimization_bounded':
    for n in range(len(Radius_beam)):
        Lattice_geom = Lattice_geometry_corrected(Lattice, Radius_beam[n])
        HomoSolid = Homogenization_3D_Solid(Lattice_geom,Cell_number,name_model,name_Job,Radius_beam[n])
        Lattice_geom = Lattice_geometry(Lattice,Radius_beam[n])
        PourcentMod = GetPourcentMod(Radius_beam[n],Lattice_geom,TypePourcentMod)
        #Initialization RadiusFactor
        RadiusFactor = 1
        min_result = minimize_scalar(FuncOptimize,bounds=(0.1, 5),method='bounded',args=(HomoSolid,Lattice_geom,name_model,name_Job,Radius_beam[n],name_Part,PourcentMod,VectorOrientation,CorrectionExteriorBeam,material_type))
        f.writelines([np.array2string(min_result.x, separator=',')+','+np.array2string(min_result.fun, separator=',')+"\n"])


f.close()
if deletevar == 1:
    for var in list(globals()):
        del globals()[var]
