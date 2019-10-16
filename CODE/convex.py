# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 13:42:22 2017

@author: Arun
"""

import os
os.chdir('your directory')


# Jarvis March O(nh) - Tom Switzer <thomas.switzer@gmail.com>

TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)


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
