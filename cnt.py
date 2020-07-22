#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 19:52:24 2020

@author: mhasan
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from math import *
class CNT:
    """ An object containing information abount the dimension of CNT inside of a Metal matrix"""
    def __init__(self,m,n,a=1.42):
        '''zigzag length:2acos30,armchir length: 2a,whrer a=1.42'''
        self.bond_length=a
        self.zigzag_length=2*a*cos(pi/6)
        self.armchair_length=2*a
        self.m=m
        self.n=n
        d=self.zigzag_length/pi*sqrt(m**2+n**2+m*n)
        cos_theta= sqrt(3)*(m+n)/2/sqrt(m**2+n**2+m*n)
        self.dia=d
        self.chirality=acos(cos_theta)*180/pi
        
    def metal_matrix_dim(self,alat,matlen):
         """determin the lengh of the CNT according to the metal matrix, c-c single bond 1.54 A, double bond 1.34 A, and triple bond 1.20 but in CNT it is 1.42 A""" 
         self.alat=alat
         self.matlen=matlen
         self.lCNT= 1.42
         n= (self.matlen*alat/2+alat/2-self.lCNT*sin(pi/3))/self.lCNT/sin(pi/3)
         self.length=round(n)

    def volume_fraction(self,n,dim2,dim3):
        """Calculate volume fraction of the CNT given the number of CNT and lateral dimension of the matrix"""
        h=3.4 # interleyar spacing of graphite
        d=self.dia
        l=self.length*self.zigzag_length/2
        L=self.alat/2*self.matlen
        v=self.matlen*dim2*dim3*(self.alat/2)**3
        v_f=(n*pi*d*h*l)/(v-n*pi*d**2*L/4)
        self.vfrac=v_f

if __name__=="__main__":
    
    #lCNT=1.42
    aAl=4.046
    #huCNT=lCNT*np.sin(np.pi/3)
    huAl=aAl/2
    m=50
    n=83
    f=lambda lCNT:n*lCNT*np.sin(np.pi/3)-m*aAl/2-huAl-lCNT*np.sin(np.pi/3)
    lCNT=fsolve(f,0)
