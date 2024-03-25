import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

from Lattice_description import *

#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Functions

#*******************************************************************************************************************
#*******************************************************************************************************************



def openFile(Lattice,Type,Option):
    if 1 in Option:
        file_path = "D:/Fichiers/70_Projet_1_Homogeneisation_Abaqus/abaqus_homogenization/Dossier_Resultats/"+Type_lattice(Lattice)+"_"+Type+"PlasticityNL.txt"
    else:
        file_path = "D:/Fichiers/70_Projet_1_Homogeneisation_Abaqus/abaqus_homogenization/Dossier_Resultats/"+Type_lattice(Lattice)+"_"+Type+"Plasticity.txt"
    with open(file_path, 'r') as file:
        data = []
        for line in file:
            line = np.fromstring(line.strip('[]\n'), sep=',')
            data.append(line)
    dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime,dataRadius = processData(data)
    return dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime,dataRadius

def openFileSolid(Lattice,Type,Option):
    if 1 in Option:
        file_path = "D:/Fichiers/70_Projet_1_Homogeneisation_Abaqus/abaqus_homogenization/Dossier_Resultats/"+Type_lattice(Lattice)+"_"+Type+"PlasticityNL.txt"
    else:
        file_path = "D:/Fichiers/70_Projet_1_Homogeneisation_Abaqus/abaqus_homogenization/Dossier_Resultats/"+Type_lattice(Lattice)+"_"+Type+"Plasticity.txt"
    with open(file_path, 'r') as file:
        data = []
        for line in file:
            line = np.fromstring(line.strip('[]\n'), sep=',')
            data.append(line)
    dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime,dataRadius = processDataSolid(data)
    return dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime,dataRadius
def processDataSolid(data):
    dataRF1 = []
    dataRF2 = []
    dataRF3 = []
    dataTime = []
    CPUTime = []
    dataU1 = []
    dataU2 = []
    dataU3 = []
    dataRadius = []
    dataRF1temp = []
    dataRF2temp = []
    dataRF3temp = []
    dataTimetemp = []
    dataU1temp = []
    dataU2temp = []
    dataU3temp = []
    dataRadiustemp = []
    CPUTimetemp = []
    idx = 0
    for i in range(len(data)):
        if idx == 0:
            dataRadiustemp.append(data[i])
        elif idx == 1:
            CPUTimetemp.append(data[i])
        elif idx==2:
            dataTimetemp.append(data[i])
        elif idx==3:
            dataRF1temp.append(data[i])
        elif idx==4:
            dataRF2temp.append(data[i])
        elif idx==5:
            dataRF3temp.append(data[i])
        elif idx==6:
            dataU1temp.append(data[i])
        elif idx==7:
            dataU2temp.append(data[i])
        elif idx==8:
            dataU3temp.append(data[i])
        if idx == 0 or idx == 1:
            idx+= 1
        elif i == len(data)-1:
            idx += 1
        else:
            try:
                if np.all(data[i]==0) and np.all(data[i-1]==0) and len(data[i])<len(data[i-1]):
                    idx += 1
                elif np.all(data[i]==0) and np.all(data[i+1]==0) and len(data[i])==len(data[i+1]):
                    idx += 1
                elif (data[i+1][0] == 0 and data[i][0]!=0) or (dataU3temp and len(data[i+1]) == 1 and len(data[i+2]) == 1 and len(data[i+3]) > 1):
                    idx+= 1
            except:
                continue
        if idx == 9:
            idx = 0
            dataRadius.append(np.ravel(dataRadiustemp))
            CPUTime.append(CPUTimetemp)
            dataTime.append(dataTimetemp)
            dataRF1.append(dataRF1temp)
            dataRF2.append(dataRF2temp)
            dataRF3.append(dataRF3temp)
            dataU1.append(dataU1temp)
            dataU2.append(dataU2temp)
            dataU3.append(dataU3temp)
            dataRF1temp = []
            dataRF2temp = []
            dataRF3temp = []
            dataTimetemp = []
            dataU1temp = []
            dataU2temp = []
            dataU3temp = []
            dataRadiustemp = []
            CPUTimetemp = []
    dataRadius = np.array(dataRadius).flatten()
    CPUTime = np.array(CPUTime).flatten()
    dataRF1 = arrayConditionnement(dataRF1)
    dataRF2 = arrayConditionnement(dataRF2)
    dataRF3 = arrayConditionnement(dataRF3)
    dataU1 = arrayConditionnement(dataU1)
    dataU2 = arrayConditionnement(dataU2)
    dataU3 = arrayConditionnement(dataU3)
    dataTime = arrayConditionnement(dataTime)
    return dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime,dataRadius

def processData(data):
    dataRF1 = []
    dataRF2 = []
    dataRF3 = []
    dataTime = []
    CPUTime = []
    dataU1 = []
    dataU2 = []
    dataU3 = []
    dataRadius = []
    dataRF1temp = []
    dataRF2temp = []
    dataRF3temp = []
    dataTimetemp = []
    dataU1temp = []
    dataU2temp = []
    dataU3temp = []
    dataRadiustemp = []
    CPUTimetemp = []
    idx = 0
    for i in range(len(data)):
        if idx == 0:
            dataRadiustemp.append(data[i])
        elif idx == 1:
            CPUTimetemp.append(data[i])
        elif idx==2:
            dataTimetemp.append(data[i])
        elif idx==3:
            dataRF1temp.append(data[i])
        elif idx==4:
            dataRF2temp.append(data[i])
        elif idx==5:
            dataRF3temp.append(data[i])
        elif idx==6:
            dataU1temp.append(data[i])
        elif idx==7:
            dataU2temp.append(data[i])
        elif idx==8:
            dataU3temp.append(data[i])
        if idx == 0 or idx == 1:
            idx+= 1
        elif i == len(data)-1:
            idx += 1
        else:
            try:
                if np.all(data[i]==0) and np.all(data[i-1]==0) and len(data[i])<len(data[i-1]):
                    idx += 1
                elif np.all(data[i]==0) and np.all(data[i+1]==0) and len(data[i])==len(data[i+1]):
                    idx += 1
                elif (data[i+1][0] == 0) or (dataU3temp and len(data[i+1]) == 1 and len(data[i+2]) == 1 and len(data[i+3]) > 1):
                    idx+= 1
            except:
                continue
        if idx == 9:
            idx = 0
            dataRadius.append(np.ravel(dataRadiustemp))
            CPUTime.append(CPUTimetemp)
            dataTime.append(dataTimetemp)
            dataRF1.append(dataRF1temp)
            dataRF2.append(dataRF2temp)
            dataRF3.append(dataRF3temp)
            dataU1.append(dataU1temp)
            dataU2.append(dataU2temp)
            dataU3.append(dataU3temp)
            dataRF1temp = []
            dataRF2temp = []
            dataRF3temp = []
            dataTimetemp = []
            dataU1temp = []
            dataU2temp = []
            dataU3temp = []
            dataRadiustemp = []
            CPUTimetemp = []
    dataRadius = np.array(dataRadius).flatten()
    CPUTime = np.array(CPUTime).flatten()
    dataRF1 = arrayConditionnement(dataRF1)
    dataRF2 = arrayConditionnement(dataRF2)
    dataRF3 = arrayConditionnement(dataRF3)
    dataU1 = arrayConditionnement(dataU1)
    dataU2 = arrayConditionnement(dataU2)
    dataU3 = arrayConditionnement(dataU3)
    dataTime = arrayConditionnement(dataTime)
    return dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime,dataRadius

def arrayConditionnement(data):
    data_combined = [np.concatenate(sublist) for sublist in data]
    return data_combined



def plotDataByRadius(dataRF_beam, dataU_beam, dataRF_beamModVar, dataU_beamModVar, dataRF_solid, dataU_solid, dataRadius, dataVolume):
    dataU_beam = [2 * arr for arr in dataU_beam]
    dataU_solid = [2 * arr for arr in dataU_solid]
    dataU_beamModVar = [2 * arr for arr in dataU_beamModVar]
    reduction_data = 1
    unique_radii = np.unique(dataRadius)
    num_radii = len(unique_radii) - reduction_data

    ncols = 3
    nrows = num_radii // ncols + (num_radii % ncols > 0)

    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 5))
    axes = axes.flatten()

    for ax in axes[num_radii:]:
        ax.remove()
    axes = axes[:num_radii]

    for i, radius in enumerate(unique_radii[:-reduction_data]):
        # Filtrer les données pour un déplacement <= 0.3 pour tous les types
        beam_indices = dataU_beam[i] <= 0.31
        beamModVar_indices = dataU_beamModVar[i] <= 0.31
        solid_indices = dataU_solid[i] <= 0.31

        # Sélectionner une valeur sur deux pour les données Solid
        solid_selection_indices = np.arange(len(dataU_solid[i])) % 2 == 0
        solid_combined_indices = np.logical_and(solid_indices, solid_selection_indices)

        filtered_dataU_beam = dataU_beam[i][beam_indices]
        filtered_dataRF_beam = dataRF_beam[i][beam_indices]

        filtered_dataU_beamModVar = dataU_beamModVar[i][beamModVar_indices]
        filtered_dataRF_beamModVar = dataRF_beamModVar[i][beamModVar_indices]

        filtered_dataU_solid = dataU_solid[i][solid_combined_indices]
        filtered_dataRF_solid = dataRF_solid[i][solid_combined_indices]

        axes[i].plot(filtered_dataU_beam, filtered_dataRF_beam, '--x', label='Beam methodology')
        axes[i].plot(filtered_dataU_beamModVar, filtered_dataRF_beamModVar, '-x', label='Our modified beam methodology')
        axes[i].plot(filtered_dataU_solid, filtered_dataRF_solid, '-o', label='Solid methodology', markersize=5)

        axes[i].set_title('Relative density : ' + str(round(dataVolume[i], 2)),fontsize = 15)
        axes[i].set_xlabel('Macroscopic strain',fontsize = 15)
        # axes[i].set_ylabel('Macroscopic stress (Pa)')
        # axes[i].legend(loc='upper left')
        axes[i].grid(True)
    axes[0].set_ylabel('Macroscopic stress (Pa)',fontsize = 15)
    axes[3].set_ylabel('Macroscopic stress (Pa)',fontsize = 15)
    plt.tight_layout()

def processDataStressStrain(dataRF, dataU,Number_cell):
    dataRF = [-rf / (Number_cell * Number_cell) for rf in dataRF]
    dataU = [-u / Number_cell for u in dataU]
    return dataRF, dataU


#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Variables

#*******************************************************************************************************************
#*******************************************************************************************************************
Case = 1
Number_cell = 1
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
Option = [1]
# Option = 0 => Nothing
# Option = 1 => Non-Linear behavior
Langue = "FR"

#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Core

#*******************************************************************************************************************
#*******************************************************************************************************************

# Beam Data
Type = 'Beam'
dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime,dataRadius = openFile(Lattice,Type, Option)
dataRF_beam, dataU_beam = processDataStressStrain(dataRF3, dataU3,Number_cell)
# BeamMod Data
Type = 'BeamModVariable'
dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime,dataRadius = openFile(Lattice,Type,Option)
dataRF_beamModVar, dataU_beamModVar = processDataStressStrain(dataRF3, dataU3,Number_cell)
# Solid Data
Type = 'Solid'
dataRF1,dataRF2,dataRF3,dataU1,dataU2,dataU3,dataTime,CPUTime,dataVolume = openFileSolid(Lattice,Type,Option)
dataRF_solid, dataU_solid = processDataStressStrain(dataRF3, dataU3,Number_cell)
plotDataByRadius(dataRF_beam, dataU_beam, dataRF_beamModVar, dataU_beamModVar, dataRF_solid, dataU_solid, dataRadius, dataVolume)

plt.show()