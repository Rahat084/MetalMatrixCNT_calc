#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 19:52:24 2020

@author: mhasan
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from math import acos
class CNT:
    """ An object containing information abount the dimension of CNT inside of a Metal matrix"""
    def __init__(self,a,m,n):
        self.bond_lengh=a
        self.m=m
        self.n=n
        d=a/np.pi*np.sqrt(m**2+n**2+m*n)
        cos_theta= np.sqrt(3)*(m+n)/2/np.sqrt(m**2+n**2+m*n)
        self.dia=d
        self.chirality=acos(cos_theta)*180/np.pi
        
    def metal_matrix_dim(self,alat,matlen):
         """determin the lengh of the CNT according to the metal matrix""" 
         self.alat=alat
         self.matlen=matlen
         self.lCNT= 1.4528
         f=lambda n: n*self.lCNT*np.sin(np.pi/3)-self.matlen*self.alat/2-self.alat/2-self.lCNT*np.sin(np.pi/3)
         self.length=round(fsolve(f,0)[0])



if __name__=="__main__":
    
    #lCNT=1.42
    aAl=4.046
    #huCNT=lCNT*np.sin(np.pi/3)
    huAl=aAl/2
    m=50
    n=83
    f=lambda lCNT:n*lCNT*np.sin(np.pi/3)-m*aAl/2-huAl-lCNT*np.sin(np.pi/3)
    lCNT=fsolve(f,0)
