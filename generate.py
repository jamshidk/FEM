#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 19:51:49 2024

@author: jamkavosi
"""

import numpy as np
import matplotlib.pyplot as plt
from uniform_mesh import *

d1=1
d2=1
p=5
m=3
R=0.2
element_type= 'D2TR3N'

#NL, EL= uniform_mesh(d1,d2,p,m, element_type)

NL, EL= void_mesh(d1,d2,p,m, R, element_type)

NoN=np.size(NL, 0)
NoE=np.size(EL, 0)

plt.figure(1)

count=1

for i in range(0, NoN):
    plt.annotate(count, xy=(NL[i,0],NL[i,1]))
    count +=1
    
    
if element_type == 'D2QU4N':
    count2=1
    
    for j in range(0,NoE):
        
        a=(NL[EL[j,0]-1,0]+NL[EL[j,1]-1,0]+NL[EL[j,2]-1,0]+NL[EL[j,3]-1,0])/4
        b=(NL[EL[j,0]-1,1]+NL[EL[j,1]-1,1]+NL[EL[j,2]-1,1]+NL[EL[j,3]-1,1])/4
        plt.annotate(count2, xy=(a,b), c='blue')
        count2 +=1
        
   
    x0,y0=NL[EL[:,0]-1,0] , NL[EL[:,0]-1,1]
    x1,y1=NL[EL[:,1]-1,0] , NL[EL[:,1]-1,1]
    x2,y2=NL[EL[:,2]-1,0] , NL[EL[:,2]-1,1]
    x3,y3=NL[EL[:,3]-1,0] , NL[EL[:,3]-1,1]

    plt.plot(np.array([x0,x1]),np.array([y0,y1]),'red',linewidth=3)
    plt.plot(np.array([x1,x2]),np.array([y1,y2]),'red',linewidth=3)
    plt.plot(np.array([x2,x3]),np.array([y2,y3]),'red',linewidth=3)
    plt.plot(np.array([x3,x0]),np.array([y3,y0]),'red',linewidth=3)
    
    plt.show()
    
if element_type == 'D2TR3N':
    count2=1
    
    for j in range(0,NoE):
        
        a=(NL[EL[j,0]-1,0]+NL[EL[j,1]-1,0]+NL[EL[j,2]-1,0])/3
        b=(NL[EL[j,0]-1,1]+NL[EL[j,1]-1,1]+NL[EL[j,2]-1,1])/3
        plt.annotate(count2, xy=(a,b), c='blue')
        count2 +=1
        
   
    x0,y0=NL[EL[:,0]-1,0] , NL[EL[:,0]-1,1]
    x1,y1=NL[EL[:,1]-1,0] , NL[EL[:,1]-1,1]
    x2,y2=NL[EL[:,2]-1,0] , NL[EL[:,2]-1,1]
    

    plt.plot(np.array([x0,x1]),np.array([y0,y1]),'red',linewidth=3)
    plt.plot(np.array([x1,x2]),np.array([y1,y2]),'red',linewidth=3)
    plt.plot(np.array([x2,x0]),np.array([y2,y0]),'red',linewidth=3)
    
    
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
