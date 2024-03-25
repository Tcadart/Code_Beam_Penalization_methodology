import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

from Lattice_description import *
#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Functions

#*******************************************************************************************************************
#*******************************************************************************************************************

def openFile(Lattice,Type):
    file_path = "D:/Fichiers/70_Projet_1_Homogeneisation_Abaqus/abaqus_homogenization/Dossier_Resultats/"+Type_lattice(Lattice)+"_"+Type+"Freq"+".txt"
    with open(file_path, 'r') as file:
        data = []
        for line in file:
            if not line.strip():
                continue
            if ',' not in line:
                value = float(line.strip())
                data.append(value)
            else:
                values = list(map(float, line.strip().split(',')))
                data.append(tuple(values))
    return data


def radiusNumber(data):
    radiusCount = 0
    radius = []
    for i in range(len(data)):
        if type(data[i]) is not tuple:
            radiusCount += 1
            radius.append(data[i])
    return radiusCount, radius

def separateData(data):
    zones = []
    taille_zone = 21

    for i in range(0, len(data), taille_zone):
        zone = data[i+1:i + taille_zone]
        zones.append(zone)
    return zones

def calculateRelativeError_norm(dataSolid,data,radius):
    relativeerror = []
    for rad in range(0, len(dataSolid)):
        vect1 = []
        vect2 = []
        for i in range(20):
            if dataSolid[rad][i][2]>1 and data[rad][i][2]>1:
                vect1.append(1/dataSolid[rad][i][2]-1/data[rad][i][2])
            else:
                vect1.append(0)
            if dataSolid[rad][i][2]>1:
                vect2.append(1/dataSolid[rad][i][2])
            else:
                vect2.append(0)
        error = LA.norm(vect1)/LA.norm(vect2)
        relativeerror.append(error)
    return relativeerror

def plotRelativeErrorFonctionDensity(dataSolide, dataBeam, dataBeamMod, dataBeamModDens,Lattice,Rad):
    radiusCount, radius = radiusNumber(dataSolide)

    dataSolideTemp = separateData(dataSolide)
    dataBeamTemp = separateData(dataBeam)
    dataBeamModTemp = separateData(dataBeamMod)
    dataBeamModDensTemp = separateData(dataBeamModDens)

    RelativeErrorBeam = calculateRelativeError_norm(dataSolideTemp, dataBeamTemp,radius)
    RelativeErrorBeamMod = calculateRelativeError_norm(dataSolideTemp, dataBeamModTemp,radius)
    RelativeErrorBeamModDens = calculateRelativeError_norm(dataSolideTemp, dataBeamModDensTemp,radius)

    RelativeErrorBeam,Rad_beam = GetLimits(Lattice,RelativeErrorBeam,radius,Rad)
    RelativeErrorBeamMod, Rad_beamMod = GetLimits(Lattice, RelativeErrorBeamMod, radius, Rad)
    RelativeErrorBeamModDens, Rad_beamModDens = GetLimits(Lattice, RelativeErrorBeamModDens, radius, Rad)

    fig = plt.plot()
    plt.semilogx(Rad_beam[1:],RelativeErrorBeam[1:],"v-",label="Model A: Beam model",linewidth = 2, markersize = 12)
    plt.semilogx(Rad_beamMod[1:], RelativeErrorBeamMod[1:], "o-", label="Model B: Beam model with strategy 4",linewidth = 2, markersize = 12)
    plt.semilogx(Rad_beamModDens[1:], RelativeErrorBeamModDens[1:], "x-", label="Model C: Beam model with mass correction",linewidth = 2, markersize = 12)
    # plt.xlabel("Relative density", fontsize=14)
    # plt.ylabel("Relative error", fontsize=14)
    plt.xticks([0.004, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7], [0.004, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7], fontsize=30)
    plt.yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8], fontsize=30)
    plt.legend(fontsize=20)
    plt.show()

def GetLimits(Lattice,RelativeError,vol_1,Rad):
    Min = 1000
    Lattice_geom = Lattice_geometry(Lattice[0], 0.1)
    for i in range(len(Lattice_geom)):
        u = [Lattice_geom[i][3]-Lattice_geom[i][0],Lattice_geom[i][4]-Lattice_geom[i][1],Lattice_geom[i][5]-Lattice_geom[i][2]]
        l = np.linalg.norm(u)
        if l<Min:
            Min = l
    Radiusmax = Min/4
    dim = len(RelativeError)
    Rad = [float(r) for r in Rad]
    for j in range(dim):
        if Rad[j]>Radiusmax:
            RelativeError[j] = np.nan
            vol_1[j] = np.nan
    return RelativeError,vol_1

#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Variables

#*******************************************************************************************************************
#*******************************************************************************************************************
Radius_beam = np.array([0.015,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2])
Lattice = [0] # 0,4,11
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
#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Core

#*******************************************************************************************************************
#*******************************************************************************************************************

for i in range(len(Lattice)):
    SolidData = openFile(Lattice[i],'Solid')
    BeamData = openFile(Lattice[i],'Beam')
    BeamModData = openFile(Lattice[i],'BeamModVariable')
    BeamModDensData = openFile(Lattice[i],'FrequencyModDens_')
    plotRelativeErrorFonctionDensity(SolidData, BeamData, BeamModData,BeamModDensData,Lattice,Radius_beam)
