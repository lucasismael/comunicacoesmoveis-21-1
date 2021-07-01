#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 02:27:27 2021

@author: lucasismael
"""

import warnings
warnings.filterwarnings("ignore")
import numpy as np
import cmath
import matplotlib.pylab as plt
from itertools import cycle
import scipy.stats as st
import os
import argparse
import yaml
import matplotlib
from random import randint   
import itertools  
    
#vtFc = [800, 900, 1800, 1900, 2100];
def findradius(dFc, dR):
    dPasso = np.ceil(dR/100);
    dRMin = dPasso;
    dIntersiteDistance = 2*np.sqrt(3/4)*dR;
    dDimX = 5*dR;
    dDimY = 6*np.sqrt(3/4)*dR;
    dPtdBm = 57;
    dPtLinear = 10**(dPtdBm/10)*1e-3;
    dSensitivity = -104;
    dHMob = 5;
    dHBs = 30;
    dAhm = 3.2*(np.log10(11.75*dHMob))**(2) - 4.97; # Fator de correção
    
    vtBs = [0];
    dOffset = np.pi/6;
    for iBs in range(2,8):
        vtBs = np.append(vtBs, dR*np.sqrt(3)*np.exp(complex(0,(iBs-2)*np.pi/3 + dOffset)));
    vtBs = vtBs + complex(dDimX/2,dDimY/2);
    
    dDimY = np.ceil(dDimY + np.mod(dDimY, dPasso));
    dDimX = np.ceil(dDimX + np.mod(dDimX, dPasso));
    
    [mtPosx, mtPosy] = np.meshgrid(np.arange(0,dDimX+1,dPasso), np.arange(0,dDimY+1,dPasso));
    
    #for dFc in vtFc:
    mtPosEachBS = np.empty([mtPosx.shape[0],mtPosx.shape[1],len(vtBs)], dtype=complex);
    mtPowerEachBSdBm = np.empty([mtPosx.shape[0],mtPosx.shape[1],len(vtBs)]);
    mtPowerFinaldBm = -np.inf*np.ones([mtPosy.shape[0],mtPosy.shape[1]]);

    for x in range(mtPosx.shape[0]):
        for y in range (mtPosx.shape[1]):
            mtPosEachBS[x,y,:] = complex(mtPosx[x,y],mtPosy[x,y]);
        
    for iBsD in range(len(vtBs)):
        mtPosEachBS[:,:,iBsD] = mtPosEachBS[:,:,iBsD]-vtBs[iBsD];
        mtDistEachBS = np.abs(mtPosEachBS[:,:,iBsD]);
        mtDistEachBS[mtDistEachBS < dRMin] = dRMin;
        #mtPldB = 69.55 + 26.16*np.log10(dFc) + (44.9 - 6.55*np.log10(dHBs))*np.log10(mtDistEachBS/1e3) - 13.82*np.log10(dHBs) - dAhm;
        mtPldB = 46.3 + 33.9*np.log10(dFc) - 13.82*np.log10(dHBs) - dAhm + (44.9 - 6.55*np.log10(dHBs))*np.log10(mtDistEachBS/1e3) + 3;
        mtPowerEachBSdBm[:,:,iBsD] = dPtdBm - mtPldB;
        mtPowerFinaldBm = np.maximum(mtPowerFinaldBm, mtPowerEachBSdBm[:,:,iBsD]);
        
    belowsensitivity = mtPowerFinaldBm[mtPowerFinaldBm < dSensitivity];
    allelements = mtPowerFinaldBm.flatten();
    dOutRate = 100*len(belowsensitivity)/len(allelements);
    return dOutRate;
    #print("Frequência da portadora = " + str(dFc) + " MHz");
    #print("Taxa de outage = " + str(dOutRate) + " %");
    #print("--------------------------------------")
    

vtFc = [800, 900, 1800, 1900, 2100];
for dFc in vtFc:
    dR = 11e3;
    outage = 10;
    dOutRate = np.inf;
    while outage <= dOutRate: 
        dOutRate = findradius(dFc, dR);
        dR = dR - 1e2;
    print("Frequência: " + str(dFc) + " MHz");
    print("Raio mínimo: " + str(dR/1000) + " km");
    print("Outage: " + str(dOutRate) + " %");
    print("-------------------------------------------");






    
