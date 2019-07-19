# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:09:38 2017

@author: Arun
"""

import os
os.chdir('C:\\Users\\Arun\\Google Drive\\Ph.D\\RAMP\\CODE')
import numpy as np
import json
#from stl import mesh
import trimesh
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import OrderedDict
with open('Input.json', 'r') as json_data:
     data = json.load(json_data)
   
BB_L= data['Input_Parameters']['Build_Size']['L']
BB_W= data['Input_Parameters']['Build_Size']['W']
BB_H= data['Input_Parameters']['Build_Size']['H']
D= data['Input_Parameters']['D']
z_t = data['Input_Parameters']['z_t']
theta = data['Input_Parameters']['theta']
RW_b = data['Input_Parameters']['RW_b']
RW_s = data['Input_Parameters']['RW_s']
HV= data['Input_Parameters']['HV']
Fill_AG = data['Input_Parameters']['Fill_AG']
AG_s = data['Input_Parameters']['AG_s']
TV_z = data['Input_Parameters']['TV_z']
T_wipe = data['Input_Parameters']['T_wipe']
Cost_b = data['Input_Parameters']['Cost_b']
Cost_s = data['Input_Parameters']['Cost_s']
T_warmup = data['Input_Parameters']['T_warmup']
T_cooldown = data['Input_Parameters']['T_cooldown']
PR = data['Input_Parameters']['PR']
C_per_KW = data['Input_Parameters']['C_per_KW']
H_per_hr = data['Input_Parameters']['H_per_hr']
TL_initial = data['Input_Parameters']['TL_initial']
TW_initial = data['Input_Parameters']['TW_initial']
Vol_b_initial = data['Input_Parameters']['Vol_b_initial']
Vol_s_initial = data['Input_Parameters']['Vol_s_initial']
CO_2_EF = data['Input_Parameters']['CO_2_EF']
CO_2_FF_b = data['Input_Parameters']['CO_2_FF_b']
CO_2_FF_s = data['Input_Parameters']['CO_2_FF_s']
CO_2_FFR_b = data['Input_Parameters']['CO_2_FFR_b']
CO_2_FFR_s = data['Input_Parameters']['CO_2_FFR_s']
Scale = data['Input_Parameters']['Scale']
m = data['Input_Parameters']['m']
Den_b = data['Input_Parameters']['Den_b']
Den_s = data['Input_Parameters']['Den_s']
Cost_misc= data['Input_Parameters']['Cost_misc']

#Read Mesh
my_mesh = trimesh.load_mesh('Turb.stl')
vert = my_mesh.vertices * (Scale)
#center to origin
gg_0=abs(min(vert[:,0]))
gg_1=abs(min(vert[:,1]))
gg_2=abs(min(vert[:,2]))
vert[:,0]=vert[:,0] + gg_0
vert[:,1]=vert[:,1] + gg_1
vert[:,2]=vert[:,2] + gg_2

#Plot data
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x =vert[:,0]
y =vert[:,1]
z =vert[:,2]
ax.scatter(x, y, z, c='r', marker='o')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()

#Layering
H = max(vert[:,2])
Layer=z_t
f=[]
P_AB=[]

def cwangle_distance(point):
        
        vect = [point[0]-origin[0], point[1]-origin[1]]
        
        lenvect = math.hypot(vect[0], vect[1])
        
        if lenvect == 0:
            return -math.pi, 0
        
        normalize = [vect[0]/lenvect, vect[1]/lenvect]
        dot_prod  = normalize[0]*refvec[0] + normalize[1]*refvec[1]     
        diff_prod = refvec[1]*normalize[0] - refvec[0]*normalize[1]     
        Ang = math.atan2(diff_prod, dot_prod)
        
        if Ang < 0:
            return 2*math.pi+Ang, lenvect
        
        return Ang, lenvect
        
def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1))) 

z_t=0

while (H/z_t >= 1):
    f=[] 
    for i in range(len(vert)):
        if ((vert[i,2]<=(z_t+0.25)) and (vert[i,2]>= (z_t-0.25))):
            f.append(vert[i,(0,1)])
            i=i+1
            
    f=np.array(f)
    
    if len(f)!=0:
            pts = f
            origin = [2, 3]
            refvec = [0, 1]
            ff=sorted(pts, key=cwangle_distance)
            ff=np.array(ff)
            x = ff[:,0]
            y = ff[:,1]
            PA = PolyArea(x,y)
            P_AB.append(PA)  
            z_t = z_t + Layer 
    else:
        P_AB.append(PA)
        z_t=z_t + Layer
        

z_t = data['Input_Parameters']['z_t']
P_ABS=[]
a=range(len(P_AB))
for i in reversed(range(len(P_AB))):
    if P_AB[i]>P_AB[i-1]:
        P_ABS.append((P_AB[i]-P_AB[i-1])+P_AB[i])


ACS_Poly = np.array(P_AB) 
ACS_sup_poly = np.array(P_ABS) 

### Transformation Equations ###

## Initial Area of Building Base $(mm^2)$ :

A_BB_ini = BB_L * BB_W

# Number of layers :

N = int(H / z_t)

## Fraction of Area of Build material filled actually by build material:

d_b = RW_b/ (RW_b + Fill_AG )

## Build Time for Part material for layer i $(min)$:

T_b_i=[]
for i in range(0,N):
    a = d_b * (ACS_Poly[i]/(RW_b * HV))  
    T_b_i.append(a)
    
## Fraction of Area of Support material filled actually by support material:

d_s = RW_s/ (RW_s + AG_s)

## Build Time for Support Material for layer i $(sec)$:

try:
    T_s_i=[]
    for i in range(0,N):
        a = d_s * (ACS_sup_poly[i])/(RW_s * HV)  
        T_s_i.append(a)
        
except IndexError:
    gotdata = 'null'
    
    
leng = (N-len(T_s_i))    
if leng >0:
        for j in range(0,leng):
            T_s_i.append(0)

#Build Time for layer i $(sec)$:

try:
    T_i=[]
    for j in range(0,N):
        a =  T_b_i[j]+ T_s_i[j]
        T_i.append(a)
        
except IndexError:
    gotdata = 'null'
    
# Time to move table in z-direction $(hrs)$:

T_zmove = z_t/TV_z

# Total Build Time $(hrs)$:

T =  (sum(T_i) + (N * T_zmove) + ((N/m) * T_wipe)) /3600

#Total Time for part material building  $(hrs)$:

T_b = (sum(T_b_i) / 3600)

#Total Time for support material building $(hrs)$:

T_s = (sum(T_s_i) / 3600)

# Volume of build material used $(mm^3)$:

V_b = d_b * sum(ACS_Poly) * z_t 

#Area of Base $(mm^2)$:
    
Base=max(ACS_Poly)

# Volume of support material used $(mm^3)$:

V_s = d_s * ((sum(ACS_sup_poly) * z_t) + (10 * Base * z_t))

# Cost of material used (\$):

C_mat =  (V_b * Cost_b) + (V_s * Cost_s)

# Total Operation time (hr):

T_o = T + (T_warmup/60) + (T_cooldown/60) 

# Energy Usage (KWhr):
 
E = T_o * PR 
 
# Energy Cost (\$):
 
C_e = E * C_per_KW 
 
# Total Heat Expelled (BTU):
 
T_heat = T_o * H_per_hr 
 
# Tip Life Remaining (hrs):
 
T_tip = TL_initial - (T_b + T_s) 
 
# Tip-wipe Assembly Life Remaining (hrs):
 
T_tip_wipe = TW_initial - (T - ((sum(T_i) + (N * T_zmove))/3600))
 
# Volume of Material Remaining $(mm^3)$:
 
VR_b = Vol_b_initial - V_b 
VR_s = Vol_s_initial - V_s
 
# Percentage of Material Remaining (%):
 
PVR_b = 1 - (V_b/Vol_b_initial) 
PVR_s = 1 - (V_s/Vol_s_initial)

## Area of Building Base used $(mm^2)$:

BB_used = 1.05 * Base
 
## Area of Building Base remaining $(mm^2)$:

A_BB_rem = A_BB_ini - BB_used
 
# Carbon Emission (Kg):
 
CO_2_emission = E * CO_2_EF 
 
# Surface Roughness ($\mu$m):

def SR_70(theta):
    thet = theta * 0.0174
    SR = 70 * (z_t/(math.cos(thet)))
    return SR

def SR_90(theta):    
    SR = 112.6 * z_t 
    return SR

def SR_70_90(theta):
    SR=0.05*((90* SR_70(70)) - (70* SR_90(90)) + theta*(SR_90(90) - SR_70(70)))
    return SR
    
if (theta >=0 and theta <=70):
    SR = SR_70(theta)
elif (theta >70 and theta < 90): 
    SR = SR_70_90(theta)
else:
    SR = 112.5 * z_t 

    
# Mass = Density * Volume (KG):

m_b= Den_b * V_b
m_s= Den_s * V_s

#Process Productivity (Kg/hr)

PP=(m_b+m_s) / (T_b + T_s)

# Specific Energy Consumption (KWhr/Kg):

SEC = (PR/PP) 
 
# Carbon Footprint Production (Kg):
 
CO_2_F_b = m_b * CO_2_FF_b   
CO_2_F_s = m_s * CO_2_FF_s
 
# Carbon Footprint Recycling (Kg):
 
CO_2_FR_b = m_b * CO_2_FFR_b  
CO_2_FR_s = m_s * CO_2_FFR_s

#Cost of Product ($):

C_p = C_mat + C_e + Cost_misc
 
# Number of products that can be manufactured (Nos):

NP= min(math.floor((Vol_b_initial/V_b)),math.floor((Vol_s_initial/V_s))) 

# Number of baseplates required (Nos):
 
N_BB = math.ceil(NP/( A_BB_ini/BB_used))

# Horizontal Swell ratio:

HSR = RW_b / D

# Write to Output JSON

Output = {
  "ID": "Fused Deposition Modelling - UMP",
  "Output_Parameters": {
    "Material": {
          "V_b": V_b,
          "V_s" : V_s,
          "VR_b": VR_b,
          "VR_s": VR_s,
          "PVR_b": PVR_b,
          "PVR_s": PVR_s,
          "m_b": m_b,
          "m_s": m_s
    },
    "Sustainability": {
        "E": E,
        "CO_2_emission": CO_2_emission,
        "SEC": SEC,
        "CO_2_F_b": CO_2_F_b,
        "CO_2_F_s": CO_2_F_s,
        "CO_2_FR_b": CO_2_FR_b,
        "CO_2_FR_s": CO_2_FR_s,
        "T_heat": T_heat
    },
    "Cost": {
       "C_mat" : C_mat,
       "C_e": C_e,
       "C_p": C_p
    },
    "Time": {
       "T ": T,
       "T_o" : T_o,
       "T_tip": T_tip,
       "T_tip_wipe": T_tip_wipe
    },
    "Miscellaneous": {
        "SR": SR,
        "HSR": HSR,
        "PP": PP,
        "NP ": NP,
        "N_BB": N_BB
                    
    }
  }
}


with open('Output.json', 'w') as json_output_data:
     json.dump(Output, json_output_data, sort_keys=True)
 