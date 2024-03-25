import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

from Lattice_description import *
#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Functions

#*******************************************************************************************************************
#*******************************************************************************************************************

def processingdata(data):
    data_process = np.zeros((len(data),12))
    for i in range(len(data)):
        data_process[i:] = np.array([1/data[i,0],1/data[i,1],1/data[i,2],1/data[i,6],1/data[i,7],1/data[i,8],-data[i,3]/data[i,0],-data[i,3]/data[i,1],-data[i,4]/data[i,0],-data[i,4]/data[i,2],-data[i,5]/data[i,1],-data[i,5]/data[i,2]])
    return data_process

def calculationRelativeError(data_process_b,data_process_s):
    dim_data = len(data_process_b)
    normfrob_s =np.zeros(dim_data)
    normfrob_dif = np.zeros(dim_data)
    Erreur_relative = np.zeros(dim_data)
    for i in range(dim_data):
        mat_rigid = np.array(
            [[data_process_s[i,0],data_process_s[i,6],data_process_s[i,8], 0, 0, 0],
            [data_process_s[i,7],data_process_s[i,1],data_process_s[i,10], 0, 0, 0],
            [data_process_s[i,9],data_process_s[i,11],data_process_s[i,2], 0, 0, 0],
            [ 0, 0, 0, data_process_s[i,3], 0, 0],
            [ 0, 0, 0, 0, data_process_s[i,4], 0],
            [ 0, 0, 0, 0, 0, data_process_s[i,5]]]
        )
        mat_rigid_inv_s = LA.inv(mat_rigid)
        normfrob_s = LA.norm(mat_rigid_inv_s,'fro')
        mat_rigid = np.array(
            [[data_process_b[i,0],data_process_b[i,6],data_process_b[i,8], 0, 0, 0],
            [data_process_b[i,7],data_process_b[i,1],data_process_b[i,10], 0, 0, 0],
            [data_process_b[i,9],data_process_b[i,11],data_process_b[i,2], 0, 0, 0],
            [ 0, 0, 0, data_process_b[i,3], 0, 0],
            [ 0, 0, 0, 0, data_process_b[i,4], 0],
            [ 0, 0, 0, 0, 0, data_process_b[i,5]]]
        )
        mat_rigid_inv_b = LA.inv(mat_rigid)

        #Calculation norm dif
        mat_rigid_dif = mat_rigid_inv_s - mat_rigid_inv_b
        normfrob_dif = LA.norm(mat_rigid_dif,'fro')
        #Calculation relative error
        Erreur_relative[i] = normfrob_dif/normfrob_s
     
    return Erreur_relative

def openFile(Lattice,Type,ResultFile):
    with open(ResultFile+Type_lattice(Lattice)+"_"+Type+".txt","r") as f:
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
    return data

def SeparateData(data):
    Radius = data[:,0]
    vol = data[:,1]
    Dim = data.shape
    if Dim[1] == 12:
        data = data[:,2:-1]
        PourcentMod = data[:,-1]
    else:
        data = data[:,2:]
        PourcentMod = np.zeros((Dim[0],1))
    return Radius,vol,data,PourcentMod

def DifDensity(vol_b,vol_s):
    return (vol_s-vol_b)


def GetRelativeError(CaseName,Lattice,ResultFile):
    RelativeError = np.zeros((9,len(Lattice)))
    DifDens = np.zeros((9,len(Lattice)))
    vol_2 = np.zeros((9,len(Lattice)))
    vol_1 = np.zeros((9,len(Lattice)))
    for n in range(len(Lattice)):
        # Open Data Beam
        data_2 = openFile(Lattice[n],CaseName[1],ResultFile)
        Radius_2,vol_2[:,n],data_2,PourcentMod = SeparateData(data_2)
        data_process_2 = processingdata(data_2)
        # Open Data Solid
        data_1 = openFile(Lattice[n],CaseName[0],ResultFile)
        Radius_1,vol_1[:,n],data_1,PourcentMod = SeparateData(data_1)
        data_process_1 = processingdata(data_1)

        # Calculation relative error
        RelativeError[:,n] = calculationRelativeError(data_process_2,data_process_1)

        # Caculation Density gap
        DifDens[:,n] = np.around(DifDensity(vol_2[:,n],vol_1[:,n]),2)

        # Calculation Rad
        Rad = list(Radius_1.astype(str))
    return RelativeError,DifDens,Rad,vol_1,vol_2,PourcentMod





def GetRelativeErrorCase(CaseName,Lattice,ResultFile):
    RelativeError,DifDens,Rad,vol_1,vol_2,PourcentMod = GetRelativeError(CaseName[0:2],Lattice,ResultFile)
    CaseName1 = [CaseName[0],CaseName[2]]
    RelativeError_1,DifDens_1,Rad,vol_1_1,vol_2_1,PourcentMod_1 = GetRelativeError(CaseName1,Lattice,ResultFile)
    CaseName2 = [CaseName[0],CaseName[3]]
    RelativeError_2,DifDens_2,Rad,vol_1_2,vol_2_2,PourcentMod_2 = GetRelativeError(CaseName2,Lattice,ResultFile)
    CaseName3 = [CaseName[0],CaseName[4]]
    RelativeError_3,DifDens_3,Rad,vol_1_3,vol_2_3,PourcentMod_3 = GetRelativeError(CaseName3,Lattice,ResultFile)
    CaseName4 = [CaseName[0],CaseName[5]]
    RelativeError_4,DifDens_4,Rad,vol_1_4,vol_2_4,PourcentMod_4 = GetRelativeError(CaseName4,Lattice,ResultFile)
    if len(CaseName)==7:
        RelativeError_5 = OpenFileOpti(CaseName[6],Lattice,ResultFile)
        RelativeError_5 = np.transpose(np.array([RelativeError_5]))
        RelativeError = np.concatenate((RelativeError,RelativeError_1,RelativeError_2,RelativeError_3,RelativeError_4,RelativeError_5),axis=1)
        DifDens = np.concatenate((DifDens,DifDens_1,DifDens_2,DifDens_3,DifDens_4,DifDens_4),axis=1)
        vol_1 = np.concatenate((vol_1,vol_1_1,vol_1_2,vol_1_3,vol_1_4,vol_1_4),axis=1)
        vol_2 = np.concatenate((vol_2,vol_2_1,vol_2_2,vol_2_3,vol_2_4,vol_2_4),axis=1)
        PourcentMod = np.concatenate((PourcentMod,PourcentMod_1,PourcentMod_2,PourcentMod_3,PourcentMod_4,PourcentMod_4),axis=1)
    else:
        RelativeError = np.concatenate((RelativeError,RelativeError_1,RelativeError_2,RelativeError_3,RelativeError_4),axis=1)
        DifDens = np.concatenate((DifDens,DifDens_1,DifDens_2,DifDens_3,DifDens_4),axis=1)
        vol_1 = np.concatenate((vol_1,vol_1_1,vol_1_2,vol_1_3,vol_1_4),axis=1)
        vol_2 = np.concatenate((vol_2,vol_2_1,vol_2_2,vol_2_3,vol_2_4),axis=1)
        PourcentMod = np.concatenate((PourcentMod,PourcentMod_1,PourcentMod_2,PourcentMod_3,PourcentMod_4),axis=1)
    return RelativeError,DifDens,Rad,vol_1,vol_2,PourcentMod


def OpenFileOpti(CaseName,Lattice,ResultFile):
    with open(ResultFile+Type_lattice(Lattice[0])+"_"+CaseName+".txt","r") as f:
        lines = f.readlines()
    data = np.array([list(map(float, line.strip().replace('[','').replace(']','').split(','))) for line in lines])
    data = np.insert(data[:,1], 0, 0)
    return data


def GetLimits(Lattice,RelativeError,vol_1,Rad):
    Min = 1000
    Lattice_geom = Lattice_geometry(Lattice, 0.1)
    for i in range(len(Lattice_geom)):
        u = [Lattice_geom[i][3]-Lattice_geom[i][0],Lattice_geom[i][4]-Lattice_geom[i][1],Lattice_geom[i][5]-Lattice_geom[i][2]]
        l = np.linalg.norm(u)
        if l<Min:
            Min = l
    Radiusmax = Min/4
    dim = RelativeError.shape
    Rad = [float(r) for r in Rad]
    for j in range(dim[0]): 
        if Rad[j]>Radiusmax:
            RelativeError[j,:] = np.nan
            vol_1[j,:] = np.nan
    return RelativeError,vol_1

#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Variables

#*******************************************************************************************************************
#*******************************************************************************************************************
Lattice = [0]
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
# Option = 0 => No option
# Option = 1 => Data reduction
ResultFile = "D:/Fichiers/70_Projet_1_Homogeneisation_Abaqus/abaqus_homogenization/Dossier_Resultats/"

#*******************************************************************************************************************
#*******************************************************************************************************************

                                    #Core

#*******************************************************************************************************************
#*******************************************************************************************************************

CaseName = ['Solid','Beam','BeamMod','BeamModVariable','BeamModVariable2','BeamModVariablecut','Optimization_bounded']
RelativeError,DifDens,Rad,vol_1,vol_2,PourcentMod = GetRelativeErrorCase(CaseName,Lattice,ResultFile)


valabsi = range(len(RelativeError))
if 1 in Option:
    RelativeError,vol_1 = GetLimits(Lattice[0],RelativeError,vol_1,Rad)
    

markers = ['o', 's', '^', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|', '_']
label = ['Without strategy','Strategy 1: L/10','Strategy 2: Radius','Strategy 3: $\\sqrt{3}\\times Radius$','Strategy 4: Radius / tan($\\gamma /2$)','Optimized coefficient']
for k in range(0,len(Lattice)*5):
    plt.semilogx(vol_1[1:,0], list(RelativeError[1:,k]),color=Color_lattice(k),label = label[k],linewidth = 2, marker=markers[k],markersize=12)
plt.legend(fontsize = 30)

plt.ylabel('Relative error between solid and beam model',fontsize = 20)
plt.xlabel('Relative density',fontsize = 30)
plt.title('Relative error difference between solid and beam model versus relative density',fontsize = 25)

plt.xticks([0.004,0.01,0.05,0.1,0.3,0.5,0.7],[0.004,0.01,0.05,0.1,0.3,0.5,0.7],fontsize = 30)
plt.yticks([0,0.1,0.2,0.3,0.4,0.5],fontsize = 30)
plt.show()



