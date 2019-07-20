# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 13:42:22 2017

@author: Arun
"""

import os
os.chdir('your directory')


# Jarvis March O(nh) - Tom Switzer <thomas.switzer@gmail.com>

TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)

#def turn(p, q, r):
#    """Returns -1, 0, 1 if p,q,r forms a right, straight, or left turn."""
#   return cmp((q[0] - p[0])*(r[1] - p[1]) - (r[0] - p[0])*(q[1] - p[1]), 0)
#
#
#    return ((q[0] - p[0])*(r[1] - p[1]) > (r[0] - p[0])*(q[1] - p[1]) -(q[0] - p[0])*(r[1] - p[1]) < (r[0] - p[0])*(q[1] - p[1]), 0)

def turn(p, q, r):
    a= -1
    b=0
    return ((a > b) - (a < b)) 
    
def _dist(p, q):
    """Returns the squared Euclidean distance between p and q."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    return dx * dx + dy * dy

def _next_hull_pt(points, p):
    """Returns the next point on the convex hull in CCW from p."""
    q = p
    for r in points:
        t = turn(p, q, r)
        if t == TURN_RIGHT or t == TURN_NONE and _dist(p, r) > _dist(p, q):
            q = r
    return q

def convex_hull(points):
    """Returns the points on the convex hull of points in CCW order."""
    hull = [min(points)]
    for p in hull:
        q = _next_hull_pt(points, p)
        if q != hull[0]:
            hull.append(q)
    return hull
points=cc
hull= convex_hull(cc)
    
     
     
   #       "T ": T,
#       "V_b": V_b,
#       "V_s" : V_s,
#       "C_mat" : C_mat,
#       "T_o" : T_o,
#       "E": E,
#        "C_e": C_e,
#        "T_heat": T_heat,
#        "T_tip": T_tip,
#        "T_tip_wipe": T_tip_wipe,
#        "VR_b": VR_b,
#        "VR_s": VR_s,
#        "PVR_b": PVR_b,
#        "PVR_s": PVR_s,
#        "A_BB_rem": A_BB_rem,
#        "CO_2_emission": CO_2_emission,
#        "SR": SR,
#        "SEC": SEC,
#        "CO_2_F_b": CO_2_F_b,
#        "CO_2_F_s": CO_2_F_s,
#        "CO_2_FR_b": CO_2_FR_b,
#        "CO_2_FR_s": CO_2_FR_s,
#        "NP ": NP,
#        "N_BB": N_BB,
#        "H_swell": 2000,
#        "Heat Expelled (BTU/hour)": 2250,
#        "Operating Sound Pressure at Printing (dBA)": 62,
#        "Nozzle Diameter": 0.381,
#        "Ns to deposit": 10,
#        "Layer Thickness": 3,
#        "Extrusion Velocity": 25,
#        "Head Travel Speed": 8,
#        "Energy Cost (per kwh)": 0.63,
#        "CO2 (per kwh)": 0.0703,
#        "Machine Cost": 20000,
#        "Wear and Tear": 1,
#        "Number of Units that can be produced": 1,
#        "Shrinkage Factor": 1,
#        "Cost of the product": 145,
#        "Layer Resolution (Slice Thickness in mm)": 0.33,
#        "Minimum Wall thichness stratasys guideline (to avoid brittleness = 4 X Layer resol)": 1.32,
#        "Model Scale": 2,
#        "Envelope Temperature": 1,
#        "Build Material Temperature": 2,
#        "Support Material Temperature": 2,
#        "Nozzle diameter": 3,
#        "Flow rate": 2,
#        "Horizontal Swell": 1,
#        "Support material initial layer (1% total mesh volume)": 5,
#        "Material Wastage": 33,
#        "Carbon footprint build material": 3,
#        "Carbon footprint support material": 2,
#        "Terpolymer of Methacrylic Acid, Styrene, and Butylacrylate": 4
#        
     
     
     
#    
#except (ValueError, TypeError):
#    print('empty list or invalid input')    
#        
        
        #### gggg    #####
#    
#    ff=[]
#    uni = set(f[:,0])
#    for i in uni:
#        
#        x = f[:,1]
#        y = f[:,0]
#    
#        # y should be sorted for both of these methods
#        order = y.argsort()
#        y = y[order]
#        x = x[order]
#        
#        def what_is_x_when_y_is(input, x, y):
#            return x[y.searchsorted(input, 'left')]
#        
#        def interp_x_from_y(input, x, y):
#            return np.interp(input, y, x)
#            
#        b= i
#        
#        a= (what_is_x_when_y_is(b, x, y))
#        ff.append([b,a])
#    
#    ff.append([0,0])
    #shoelace formula
#    ff=np.array(ff)
#    from operator import itemgetter
#    ff=np.array(sorted(ff, key=itemgetter(0)))
#    x = ff[:,0]
#    y = ff[:,1]
# 
#import math
#
#pts = f
#origin = [2, 3]
#refvec = [0, 1]
#
#def clockwiseangle_and_distance(point):
#    # Vector between point and the origin: v = p - o
#    vector = [point[0]-origin[0], point[1]-origin[1]]
#    # Length of vector: ||v||
#    lenvector = math.hypot(vector[0], vector[1])
#    # If length is zero there is no angle
#    if lenvector == 0:
#        return -math.pi, 0
#    # Normalize vector: v/||v||
#    normalized = [vector[0]/lenvector, vector[1]/lenvector]
#    dotprod  = normalized[0]*refvec[0] + normalized[1]*refvec[1]     # x1*x2 + y1*y2
#    diffprod = refvec[1]*normalized[0] - refvec[0]*normalized[1]     # x1*y2 - y1*x2
#    angle = math.atan2(diffprod, dotprod)
#    # Negative angles represent counter-clockwise angles so we need to subtract them 
#    # from 2*pi (360 degrees)
#    if angle < 0:
#        return 2*math.pi+angle, lenvector
#    # I return first the angle because that's the primary sorting criterium
#    # but if two vectors have the same angle then the shorter distance should come first.
#    return angle, lenvector
#
#ff=sorted(pts, key=clockwiseangle_and_distance)
#ff=np.array(ff)
#
#x = ff[:,0]
#y = ff[:,1]
#    
#def PolyArea(x,y):
#    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
#
#PA = PolyArea(x,y)
#P_AB.append(PA)  
#z_t = z_t + Layer  
#    #### gggg    #####
#    
#    
#    bead=25
#    b_width = max(f[:,1])
#    while (b_width/bead >= 1):
#        g=[]
#        for i in range(len(f)):
#            if ((f[i,1]<=bead) and (f[i,1]>= (bead-Bead))):
#                g.append(f[i,0])
#                i=i+1
#        
#        g=np.array(g)
#        max_g=max(g)
#        min_g=min(g)
#        dist = max_g - min_g
        
#        bead = bead + Bead
      
     
     
     
     
     
     
     
     

#
##f=np.sort(f, axis=0)
#
##new_array = [tuple(column) for column in f]
##uniques = np.unique(new_array)
##import operator
##ff=f[:,(0,1)]
##ff=np.around(ff, decimals=2)
##copy = []
##[
##    copy.append(max(
##        (point2 for point2 in ff if point["x"] == point2["x"]),
##        key=operator.itemgetter('y'))
##    )
##    for point in ff
##        if next(
##            (cPoint for cPoint in copy if cPoint['x'] == point['x']),
##            None) == None
##]
##maxima = {}
##for d in ff:
##    x, y = d[0], d[1]
##    if x not in maxima or y > maxima[x]:
##        maxima[x] = y
##copy = [d for d in lst if d[1] == maxima[d[0]]]
#
#
##col=0
##ff=ff[np.argsort(ff[:,col])]
#
##
##def unique_rows(f):
##    f = np.ascontiguousarray(f)
##    unique_f = np.unique(f.view([('', f.dtype)]*f.shape[1]))
##    return unique_f.view(f.dtype).reshape((unique_f.shape[0], f.shape[1]))
##
##uniques= unique_rows(ff)
##ff=uniques
#
#####   GOOD   ####
#ff=[]
#uni = set(f[:,0])
#for i in uni:
#    
#    x = f[:,1]
#    y = f[:,0]
#
#    # y should be sorted for both of these methods
#    order = y.argsort()
#    y = y[order]
#    x = x[order]
#    
#    def what_is_x_when_y_is(input, x, y):
#        return x[y.searchsorted(input, 'left')]
#    
#    def interp_x_from_y(input, x, y):
#        return np.interp(input, y, x)
#        
#    b= i
#    
#    a= (what_is_x_when_y_is(b, x, y))
#    ff.append([b,a])
#
#ff.append([0,0])
##shoelace formula
#ff=np.array(ff)
#from operator import itemgetter
#ff=np.array(sorted(ff, key=itemgetter(0)))
#x = ff[:,0]
#y = ff[:,1]
#
#def PolyArea(x,y):
#    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
#
#math.pi*50*50*10    
#print (PolyArea(x,y))
#
##### gggg    #####
#Initial_contour= 3*Layer*1.1*PolyArea(x,y)
#
#Part_Support = my_mesh.convex_hull.volume - my_mesh.volume
#Total_Support= Initial_contour + Part_Support
#
##area of building base:
#    Initial_Area = 203*203 
#    BB_used= 1.1*PolyArea(x,y)
#    Area_left= Initial_area - BB_used
#   
##ACS Nozzle (Assumed to be circular)
#    ACS_Nozzle = math.pi * (Nozzle_Dia/2) ** 2 
#   
##VFR
#   VFR= ACS_Nozzle * Extrusion_Velocity 
#   
##Build_Time
#    Build_Time = sum(Distance) / Speed_head
#    Build_support_time = Total_Support / VFR_support   
#    Miscellaneous_time = 10% 
#    Total = (Build_Time + Build_support_time)*1.1
#        
#
##Volume used 
#    V=  VFR * Build_Time
#   
## Cost of Build Material:
#    Build_Material_Cost (per cu mm)= 50
#    Cost_build = Build_Material_Cost * Volume used
##Energy usage
#    energy = (Warmup+Total_time + Cooldown) * power rating 
#    CO2 = 
#    Water=
#    Electricity cost =
#    Heat Expelled = 
#
##Tip Life
#    Tip_Life = Initial - (Build_Time + Build_support_time)
#    Tip_wipe = Initial_wipe -  (Total - (Build_Time + Build_support_time))
#
##Material Remaining
#Build_Material_Remaining = Initial - Volume_used
#Support_Material_Rem = Initial - Volume_used
#No_of_parts_mfg = Build_Material_Remaining / Volume_used
#
##Volume after shrinkage
#VAS= my_mesh.volume * (1-0.01)
#
##Carbon footprint 1.5 kg/kg
#Carbon_footprint_Build_production= volume(V)*1.5
##Water_usage_production= 185l/kg
#Carbon_footprint_Build_recycle= volume(V)*1.3
#Carbon_footprint_Support= volume *3
#
#Carbon_emission= energy* Carbon_intensity_factor
#
##Surface Roughness Micro meter
#Build_Orientation= 50
#Surface_Roughness= 70*(Layer/ math.cos(Build_Orientation))
#
##Specific energy consumption MJ/Kg
#
#SEC = Power_Rate_(KW)/ ProcessProductivity_(Kg/h)
#
#SEC*0.27777 #KWhr
#
#
#
#
#
##
##     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
##     
##     
##     
##dd = {'key':'value'}
##print (dd)
##dia=150
##dd['dia'] = dia
##print (dd)
##
##
##import itertools
##
##getx, gety = f[0],f[1] # or use operator.itemgetter
##groups = itertools.groupby(sorted(f, key=getx), key=getx)
##m = [max(b, key=gety) for a,b in groups]
##copy = [l for l in f if l in m]
##
##
##
##
##
##import itertools
##
##getx, gety = lambda a: a['x'], lambda a: a['y'] # or use operator.itemgetter
##groups = itertools.groupby(sorted(lst, key=getx), key=getx)
##m = [max(b, key=gety) for a,b in groups]
##copy = [l for l in lst if l in m]
##
##
##import numpy as np
##
##
##
def cwangle_distance(point):
        # Vector between point and the origin: v = p - o
        vect = [point[0]-origin[0], point[1]-origin[1]]
        # Length of vector: ||v||
        lenvect = math.hypot(vect[0], vect[1])
        # If length is zero there is no angle
        if lenvect == 0:
            return -math.pi, 0
        # Normalize vector: v/||v||
        normalize = [vect[0]/lenvect, vect[1]/lenvect]
        dot_prod  = normalize[0]*refvec[0] + normalize[1]*refvec[1]     # x1*x2 + y1*y2
        diff_prod = refvec[1]*normalize[0] - refvec[0]*normalize[1]     # x1*y2 - y1*x2
        Ang = math.atan2(diff_prod, dot_prod)
        # Negative angles represent counter-clockwise angles so we need to subtract them 
        # from 2*pi (360 degrees)
        if Ang < 0:
            return 2*math.pi+Ang, lenvect
        # I return first the angle because that's the primary sorting criterium
        # but if two vectors have the same angle then the shorter distance should come first.
        return Ang, lenvect
