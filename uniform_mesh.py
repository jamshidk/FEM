#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 19:59:39 2024

@author: jamkavosi
"""


import numpy as np

def uniform_mesh(d1, d2, p, m, element_type):
    
    PD=2
    q=np.array([[0,0],[d1,0],[0,d2],[d1,d2]])  # 4 corners
    
    NoN=(p+1)*(m+1)
    NoE=p*m
    NPE=4 ## if the element type is 
    
    
    NL=np.zeros([NoN, PD])
    
    
    a= (q[1,0]-q[0,0])/p ## increment in horizontal direction
    b= (q[2,1]-q[0,1])/m ## increment in vertical direction
    
    
    
    n=0 # this will allow us to go through the NL
    
    for i in range(1, m+2):
        for j in range(1, p+2):
            
            NL[n,0]=q[0,0]+(j-1)*a
            NL[n,1]=q[0,1]+(i-1)*b
            
            n += 1
            
            
    EL=np.zeros([NoE, NPE])
    
    for i in range(1, m+1):
        for j in range(1, p+1):
            
            if j ==1:    
                EL[(i-1)*p+j-1,0]=(i-1)*(p+1)+j
                EL[(i-1)*p+j-1,1]=EL[(i-1)*p+j-1,0]+1
                EL[(i-1)*p+j-1,3]=EL[(i-1)*p+j-1,0]+p+1
                EL[(i-1)*p+j-1,2]=EL[(i-1)*p+j-1,3]+1
            
            else:
                EL[(i-1)*p+j-1,0]=EL[(i-1)*p+j-2,1]
                EL[(i-1)*p+j-1,3]=EL[(i-1)*p+j-2,2]
                EL[(i-1)*p+j-1,1]=EL[(i-1)*p+j-1,0]+1
                EL[(i-1)*p+j-1,2]=EL[(i-1)*p+j-1,3]+1
                
    
    if element_type == 'D2TR3N':
        
        NPE_new=3
        NoE_new=2*NoE
        EL_new=np.zeros([NoE_new,NPE_new])
        
        for i in range(1, NoE+1):
            #for the first trianle element
            EL_new[2*(i-1),0]=EL[i-1,0]
            EL_new[2*(i-1),1]=EL[i-1,1]
            EL_new[2*(i-1),2]=EL[i-1,2]
            #for the second trianle element
            EL_new[2*(i-1)+1,0]=EL[i-1,0]
            EL_new[2*(i-1)+1,1]=EL[i-1,2]
            EL_new[2*(i-1)+1,2]=EL[i-1,3]        
            
        EL=EL_new
        
    EL=EL.astype(int)  
         
    return NL, EL