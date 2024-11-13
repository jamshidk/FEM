#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:04:06 2024

@author: jamkavosi
"""

import numpy as np
import math

def stiffness(x, GPE):
    
    NPE=np.size(x,0)
    PD=np.size(x,1)
    
    K = np.zeros([NPE*PD,NPE*PD])
    
    coor=x.T
    
    for i in range(1, NPE+1):
        for j in range(1, NPE+1):
            
            k=np.zeros([PD,PD])
            for gp in range (1, GPE+1):
                J=np.zeros([PD,PD])
                grad=np.zeros([PD,NPE])
                (xi,eta, alpha)= GaussPoint(NPE, GPE, gp)
                grad_nat= grad_N_Nat(NPE, xi, eta)
                J= coor @ grad_nat.T
                
                grad= np.linalg.inv(J).T @ grad_nat
                
                for a in range(1, PD+1):
                    for c in range(1, PD+1):
                        for b in range(1, PD+1):
                            for d in range(1, PD+1):
                                
                                k[a-1,c-1] = k[a-1,c-1]+ grad[b-1, i-1]*constitutive(a,b,c,d)*grad[d-1,j-1]*np.linalg.det(J)*alpha
                                
            K[((i-1)*PD+1)-1:i*PD,((j-1)*PD+1)-1:j*PD ]= k
    return K
                                

def GaussPoint(NPE, GPE, gp):
    if (NPE ==4):
        if (GPE==1):
            if (gp==1):
                xi=0
                eta=0
                alpha=4
                
        if (GPE==4):
            if (gp==1):
                xi= -1/math.sqrt(3)
                eta= -1/math.sqrt(3)
                alpha=1
                
            if (gp==2):
                xi=1/math.sqrt(3)
                eta=-1/math.sqrt(3)
                alpha=1
                
            if (gp==3):
                xi=1/math.sqrt(3)
                eta=1/math.sqrt(3)
                alpha=1
                
            if (gp==4):
                xi=-1/math.sqrt(3)
                eta=1/math.sqrt(3)
                alpha=1
    if (NPE==3):
        if (GPE==1):
            if (gp==1):
                xi=1/3
                eta=1/3
                alpha=1
            
    return (xi,eta, alpha)

def grad_N_Nat(NPE, xi, eta):
    
    PD=2
    result=np.zeros([PD, NPE])
    
    if (NPE==3):
        
        result[0,0]=  1
        result[0,1]=  0
        result[0,2]= -1
        
        
        result[1,0]=  0
        result[1,1]=  1
        result[1,2]= -1

    
    if (NPE==4):
        
        result[0,0]= -1/4*(1-eta)
        result[0,1]=  1/4*(1-eta)
        result[0,2]=  1/4*(1+eta)
        result[0,3]= -1/4*(1+eta)
        
        result[1,0]= - 1/4*(1-xi)
        result[1,1]=  -1/4*(1+xi)
        result[1,2]=   1/4*(1+xi)
        result[1,3]=   1/4*(1-xi)
        
    return result

def constitutive(i,j,k,l):
    E=8/3
    Nu=1/3
    C= (E/(2*(1+Nu)))*(delta(i,l)*(delta(j,k))+(delta(i,k)*delta(j,l))) +((E*Nu)/(1-Nu**2))*delta(i,j)*delta(k,l)
    
    return C

def delta(i,j):
    if (i==j):
        return 1
    else:
        return 0
    
    