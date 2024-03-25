import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize

from Lattice_description import *


#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Functions

#*******************************************************************************************************************
#*******************************************************************************************************************


def openFile_bounded(Lattice,DataFile):
    with open(DataFile+Type_lattice(Lattice)+"_"+"Optimization_bounded"+".txt","r") as f:
        file = f.readlines()
    Data = []
    for line in file:
        line = line.strip('[]\n').replace('[', '').replace(']', '')
        numbers = line.split(',')
        array = np.array(numbers, dtype=float)
        Data.append(array)
    Data = np.array(Data)
    return Data

def openFiledensity(Lattice,Type,DataFile):
    with open(DataFile+Type_lattice(Lattice)+"_"+Type+".txt","r") as f:
        file = f.readlines()
    arrays = []
    for line in file:
        array = np.fromstring(line.strip('[]\n'), sep=',')
        arrays.append(array)
    arrays = np.concatenate(arrays, axis=0)
    if len(arrays)%11==0:
        data = np.reshape(arrays,(int(len(arrays)/11),11))
    if len(arrays)%12==0:
        data = np.reshape(arrays,(int(len(arrays)/12),12))
    Radius = data[:,0]
    vol = data[:,1]
    Dim = data.shape
    if Dim[1] == 12:
        data = data[:,2:-1]
        PourcentMod = data[:,-1]
    else:
        data = data[:,2:]
        PourcentMod = np.zeros((Dim[0],1))
    return vol



def GetLimits(Lattice,RelativeDensity,Data,Radius_beam):
    Min = 1000
    Lattice_geom = Lattice_geometry(Lattice, 0.1)
    for i in range(len(Lattice_geom)):
        u = [Lattice_geom[i][3]-Lattice_geom[i][0],Lattice_geom[i][4]-Lattice_geom[i][1],Lattice_geom[i][5]-Lattice_geom[i][2]]
        l = np.linalg.norm(u)
        if l<Min:
            Min = l
    Radiusmax = Min/4
    dim = RelativeDensity.shape
    for j in range(dim[0]-1):
        if Radius_beam[j]>Radiusmax:
            RelativeDensity[j+1] = np.nan
            Data[j] = np.nan
    return RelativeDensity,Data

#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Variables

#*******************************************************************************************************************
#******************************************************************************************************************
Case = 5
#Case 1 => Only Radius increase
#Case 2 => Graphique 2 en 1
#Case 3 => Study of optimization
#Case 4 => Comparaison Optimization
#Case 5 => Optimization Brent / Value radius modified
Lattice = [0,4,9,3,6]
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

Radius_beam = np.array([0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2])
DataFile = "D:/Fichiers/70_Projet_1_Homogeneisation_Abaqus/abaqus_homogenization/Dossier_Resultats/"
#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Core

#*******************************************************************************************************************
#*******************************************************************************************************************

for i in range(len(Lattice)):
    # Open files with data
    Data = openFile_bounded(Lattice[i],DataFile)
    RelativeDensity = openFiledensity(Lattice[i],'Solid',DataFile)
    #Reduce domain of study
    RelativeDensity,Data = GetLimits(Lattice[i],RelativeDensity,Data,Radius_beam)
    label = []
    for n in range(len(Lattice)):
        label.append(Type_lattice(Lattice[n]))
    markers = ['o', 's', '^', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|', '_']
    plt.plot(RelativeDensity[1:],Data[:,0],label = label[i],linewidth=3,c=Color_lattice(Lattice[i]),marker = markers[i],markersize=12)
plt.title('Optimal radius increase versus relative density',fontsize = 25)
plt.xlabel('Relative density',fontsize = 30)
plt.ylabel('Radius increase at node',fontsize = 30)
plt.plot([0.001, 0.55], [1.5, 1.5], label='Proposed coefficient', linewidth=3, c='k')

plt.xticks(fontsize = 30)
plt.yticks(fontsize = 30)
plt.legend()
plt.show()