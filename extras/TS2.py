#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 19:55:45 2022

@author: bryn
"""
import numpy as np
from scipy.integrate import solve_bvp

# The Sasaki metric on TS2 in spherical coordinates
def metric(Y):
    G = np.zeros((4,4))
    phi, theta, vphi, vtheta = Y
    s = np.sin(phi)
    c = np.cos(phi)
    g11 = 1+(s*vtheta)**2
    
    G=np.array([[g11, s**2*vphi*vtheta, 0, -c*s*vtheta],
                [s**2*vphi*vtheta, c**2*g11+s**2*vphi**2, c*s*vtheta, -c*s*vphi],
                [0, c*s*vtheta, 1, 0],
                [-c*s*vtheta, -c*s*vphi, 0, c**2]])
    return G


# The Levi Civita connection on TS2 expressed in spherical coordinates
def Christoffel(Y):
    phi, theta, vphi, vtheta = Y

    s=np.sin(phi)
    c=np.cos(phi)
    t=np.tan(phi)

    Gamma = np.zeros((4,4,4))
    Gamma[0,0,1] = -s*c*vtheta*vphi/2
    Gamma[0,1,0] = Gamma[0,0,1]
    Gamma[0,1,1] = -c*s*(c**2*vtheta**2+vphi**2-1)
    Gamma[0,1,2] = -c**2*vtheta/2
    Gamma[0,2,1] = Gamma[0,1,2]
    Gamma[0,1,3] = c**2*vphi/2
    Gamma[0,3,1] = Gamma[0,1,3]
    Gamma[1,0,0] = t*vtheta*vphi
    Gamma[1,0,1] = t*(c**2*vtheta**2+vphi**2-2)/2
    Gamma[1,1,0] = Gamma[1,0,1]
    Gamma[1,0,2] = vtheta/2
    Gamma[1,2,0] = Gamma[1,0,2]
    Gamma[1,0,3] = -vphi/2
    Gamma[1,3,0] = Gamma[1,0,3]
    Gamma[2,0,0] = -s**2*vtheta**2*vphi
    Gamma[2,0,1] = (c**4*vtheta**2+(vphi**2-vtheta**2+1)*c**2-vphi**2)*vtheta/2
    Gamma[2,1,0] = Gamma[2,0,1]
    Gamma[2,0,2] = -c*s*vtheta**2/2
    Gamma[2,2,0] = Gamma[2,0,2]
    Gamma[2,0,3] = s*c*vtheta*vphi/2
    Gamma[2,3,0] = Gamma[2,0,3]
    Gamma[2,1,1] = -s**2*vphi
    Gamma[2,1,3] = c*s
    Gamma[2,3,1] = Gamma[2,1,3]
    Gamma[3,0,0] = -vtheta*(c**2*vphi**2+c**2-vphi**2)/c**2
    Gamma[3,0,1] = -vphi*((vphi**2-1)*c**2-vphi**2+2)/2/c**2
    Gamma[3,1,0] = Gamma[3,0,1]
    Gamma[3,0,2] = t*vtheta*vphi/2
    Gamma[3,2,0] = Gamma[3,0,2]
    Gamma[3,0,3] = -s*(vphi**2+2)/2/c
    Gamma[3,3,0] = Gamma[3,0,3]
    Gamma[3,1,1] = -s**2*vtheta*(c**2*vtheta**2+vphi**2)
    Gamma[3,1,2] = -t*(c**2*vtheta**2+2)/2
    Gamma[3,2,1] = Gamma[3,1,2]
    Gamma[3,1,3] = s*vtheta*vphi*c/2
    Gamma[3,3,1] = Gamma[3,1,3]
    
    return Gamma

# The second order system for the geodesics in spherical coordinates
def geodesiceq(Z):
    Y = Z[:4]
    V = Z[4:]
    Gamma = Christoffel(Y)
    Zd = np.zeros((8,))
    Zd[:4] = V
    Zd[4]  = -V.transpose()@Gamma[0,:,:]@V
    Zd[5]  = -V.transpose()@Gamma[1,:,:]@V
    Zd[6]  = -V.transpose()@Gamma[2,:,:]@V
    Zd[7]  = -V.transpose()@Gamma[3,:,:]@V
    
    return Zd

# This is a wrapper that adapts geodesiceq to the requirements of the scipy solver solve_bvp
def geo_func(x,Z):
    n, m = Z.shape
    Z1 = np.zeros_like(Z)
    for k in range(m):
        Z1[:,k] = geodesiceq(Z[:,k])
    return Z1

# This function calculates the geodesic that connects two points on TS2 for the Sasaki metric in spherical coords
# It also calculates the distance between the points
def geodesic_bvp(Y0,Y1,N=101):
    x = np.linspace(0., 1., N)
    ZZ = np.zeros((8,N))
    YY = np.column_stack((Y0,Y1))@np.vstack((1-x,x))
    ZZ[:4,:]=YY
    
    boundarycondition.Y0=Y0
    boundarycondition.Y1=Y1
    Result = solve_bvp(geo_func, boundarycondition, x, ZZ)
    L = 0
    m=Result.x.size-1
    s=Result.x
    Z=Result.y
    for k in range(m):
        Vm = 0.5*(Z[4:,k]+Z[4:,k+1])
        Ym = 0.5*(Z[:4,k]+Z[:4,k+1])
        L += (s[k+1]-s[k])*np.sqrt(Vm.transpose()@metric(Ym)@Vm)
    
    return Result,L

    
# A function defining the boundary conditions for the geodesics eqns. Required by solve_bvp    
def boundarycondition(Za, Zb):
    R = np.zeros((8,))
    R[:4] = Za[:4]-boundarycondition.Y0
    R[4:] = Zb[:4]-boundarycondition.Y1
    return R

# Transform the spherical coordinates on TS2 to cartesian coordinates
def sph2vecs(Y):
    x,y,z = sph2cart(Y[1], Y[0], 1)
    q = np.array([x,y,z])
    sp = np.sin(Y[0])
    cp = np.cos(Y[0])
    st = np.sin(Y[1])
    ct = np.cos(Y[1])
    v = np.array([-sp*ct*Y[2]-cp*st*Y[3], -sp*st*Y[2]+cp*ct*Y[3],cp*Y[2]])
    return q,v

# Transform spherical coordinates to cartesian coordinates on S2
def sph2cart(az, el, r):
    rcos_phi = r * np.cos(el)
    x = rcos_phi * np.cos(az)
    y = rcos_phi * np.sin(az)
    z = r * np.sin(el)
    return x, y, z
