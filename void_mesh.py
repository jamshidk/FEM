#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 19:37:06 2024

@author: jamkavosi
"""

import numpy as np
import math

def void_mesh(d1, d2, p, m, R, element_type):
    
    q=np.array([[0,0],[d1,0],[0,d2],[d1,d2]])
    
    PD=2
    NoN=2*(p+1)*(m+1)+2*(p-1)*(m+1)
    
    NoE=4*p*m
    NPE= 4 #D2QU4N
    
    ### Node ###
    
    NL=np.zeros([NoN, PD])
    
    a= (q[1,0]-q[0,0])/p ## increment in horizontal direction
    b= (q[2,1]-q[0,1])/p ## increment in vertical direction
    
    ### region 1 ###
    
    coor11=np.zeros([(p+1)*(m+1),PD])
    
    for i in range(1, p+2):
        
        coor11[i-1,0]= q[0,0]+(i-1)*a
        coor11[i-1,1]= q[0,1]
        
    for i in range(1, p+2):
        
        coor11[m*(p+1)+i-1,0]=R*np.cos((5*math.pi/4)+(i-1)*(math.pi/(2*p)))+d1/2
        coor11[m*(p+1)+i-1,1]=R*np.sin((5*math.pi/4)+(i-1)*(math.pi/(2*p)))+d2/2
        
    for i in range(1, m):
        for j in range(1, p+2):
            dx=(coor11[m*(p+1)+j-1,0]-coor11[j-1,0])/m
            dy=(coor11[m*(p+1)+j-1,1]-coor11[j-1,1])/m
            
            coor11[i*(p+1)+j-1,0]=coor11[(i-1)*(p+1)+j-1,0]+dx
            coor11[i*(p+1)+j-1,1]=coor11[(i-1)*(p+1)+j-1,1]+dy
            
    ### region 2 ###
    
    coor22=np.zeros([(p+1)*(m+1),PD])
    
    for i in range(1, p+2):
        
        coor22[i-1,0]= q[2,0]+(i-1)*a
        coor22[i-1,1]= q[2,1]
        
    for i in range(1, p+2):
        
        coor22[m*(p+1)+i-1,0]=R*np.cos((3*math.pi/4)-(i-1)*(math.pi/(2*p)))+d1/2
        coor22[m*(p+1)+i-1,1]=R*np.sin((3*math.pi/4)-(i-1)*(math.pi/(2*p)))+d2/2
    
    for i in range(1, m):
        for j in range(1, p+2):
            dx=(coor22[m*(p+1)+j-1,0]-coor22[j-1,0])/m
            dy=(coor22[m*(p+1)+j-1,1]-coor22[j-1,1])/m
            
            coor11[i*(p+1)+j-1,0]=coor22[(i-1)*(p+1)+j-1,0]+dx
            coor11[i*(p+1)+j-1,1]=coor22[(i-1)*(p+1)+j-1,1]+dy
            
    ### region 3 ###
    coor33=np.zeros([(p-1)*(m+1),PD])

    for i in range(1, p):
        coor33[i-1,0]= q[0,0]
        coor33[i-1,1]= q[0,1]+i*b
   
    for i in range(1, p):
       
       coor33[m*(p-1)+i-1,0]=R*np.cos((5*math.pi/4)-i*(math.pi/(2*p)))+d1/2
       coor33[m*(p-1)+i-1,1]=R*np.sin((5*math.pi/4)-i*(math.pi/(2*p)))+d2/2
       
       
    for i in range(1, m):
        for j in range(1, p):
            dx=(coor33[m*(p-1)+j-1,0]-coor33[j-1,0])/m
            dy=(coor33[m*(p-1)+j-1,1]-coor33[j-1,1])/m
            
            coor33[i*(p-1)+j-1,0]=coor33[(i-1)*(p-1)+j-1,0]+dx
            coor33[i*(p-1)+j-1,1]=coor33[(i-1)*(p-1)+j-1,1]+dy
            
            
    ### region 4 ###
    coor44=np.zeros([(p-1)*(m+1),PD])

    for i in range(1, p):
        coor44[i-1,0]= q[1,0]
        coor44[i-1,1]= q[1,1]+i*b
   
    for i in range(1, p):
       
       coor44[m*(p-1)+i-1,0]=R*np.cos((7*math.pi/4)+i*(math.pi/(2*p)))+d1/2
       coor44[m*(p-1)+i-1,1]=R*np.sin((7*math.pi/4)+i*(math.pi/(2*p)))+d2/2
       
       
    for i in range(1, m):
        for j in range(1, p):
            dx=(coor44[m*(p-1)+j-1,0]-coor44[j-1,0])/m
            dy=(coor44[m*(p-1)+j-1,1]-coor44[j-1,1])/m
            
            coor44[i*(p-1)+j-1,0]=coor44[(i-1)*(p-1)+j-1,0]+dx
            coor44[i*(p-1)+j-1,1]=coor44[(i-1)*(p-1)+j-1,1]+dy
            
            
    #### RE-ORDERING THE NODES ####
    
    ### ELEMENT ####
    