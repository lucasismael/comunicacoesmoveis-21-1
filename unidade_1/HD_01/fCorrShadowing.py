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

def fCorrShadowing (mtPoints, dShad, dAlphaCorr, dSigmaShad, dDimXOri, dDimYOri):
    
    dDimYS = np.ceil(dDimYOri + np.mod(dDimYOri, dShad));
    dDimXS = np.ceil(dDimXOri + np.mod(dDimXOri, dShad));
    [mtPosxshad, mtPosyshad] = np.meshgrid(np.arange(0,dDimXS+1,dShad), np.arange(0,dDimYS+1,dShad));
    mtPosShad = mtPosxshad + 1j*mtPosyshad;
    
    mtShadowingSamples = np.empty([mtPosyshad.shape[0], mtPosyshad.shape[1],8]);
    for iMap in range (0,8):
        mtShadowingSamples[:,:,iMap] = dSigmaShad*np.random.randn(mtShadowingSamples.shape[0], mtShadowingSamples.shape[1]);
    
    dSizel = mtPoints.shape[0];
    dSizec = mtPoints.shape[1];
    
    mtShadowingCorr = np.empty([dSizel,dSizec,8]);
    
    for il in range(0, dSizel):
        for ic in range(0, dSizec):
            
            dShadPoint = mtPoints[il,ic];
            
            dXIndexP1 = dShadPoint.real/dShad;
            dYIndexP1 = dShadPoint.imag/dShad;
            
            if (np.mod(dXIndexP1,1) == 0 and np.mod(dYIndexP1,1) == 0):
                dXIndexP1 = np.floor(dXIndexP1);
                dYIndexP1 = np.floor(dYIndexP1);
                dShadowingC = mtShadowingSamples[int(dYIndexP1),int(dXIndexP1),7];
                for iMap in range(0,7):
                    dShadowingERB = mtShadowingSamples[int(dYIndexP1), int(dXIndexP1), iMap];
                    mtShadowingCorr[il,ic,iMap] = np.sqrt(dAlphaCorr)*dShadowingC + np.sqrt(1-dAlphaCorr)*dShadowingERB;
            else:
                dXIndexP1 = np.floor(dXIndexP1);
                dYIndexP1 = np.floor(dYIndexP1);
                if (dXIndexP1 == mtPosyshad.shape[1] and dYIndexP1 == mtPosyshad.shape[0]):
                    dXIndexP2 = dXIndexP1-1;
                    dYIndexP2 = dYIndexP1;
                    dXIndexP4 = dXIndexP1-1;
                    dYIndexP4 = dYIndexP1-1;
                    dXIndexP3 = dXIndexP1;
                    dYIndexP3 = dYIndexP1-1;
                elif (dXIndexP1 == mtPosyshad.shape[1] - 1):
                    dXIndexP2 = dXIndexP1-1;
                    dYIndexP2 = dYIndexP1;
                    dXIndexP4 = dXIndexP1-1;
                    dYIndexP4 = dYIndexP1+1;
                    dXIndexP3 = dXIndexP1;
                    dYIndexP3 = dYIndexP1+1;
                elif (dYIndexP1 == mtPosyshad.shape[0] - 1):
                    dXIndexP2 = dXIndexP1+1;
                    dYIndexP2 = dYIndexP1;
                    dXIndexP4 = dXIndexP1+1;
                    dYIndexP4 = dYIndexP1-1;
                    dXIndexP3 = dXIndexP1;
                    dYIndexP3 = dYIndexP1-1;
                else:
                    dXIndexP2 = dXIndexP1+1;
                    dYIndexP2 = dYIndexP1;
                    dXIndexP4 = dXIndexP1+1;
                    dYIndexP4 = dYIndexP1+1;
                    dXIndexP3 = dXIndexP1;
                    dYIndexP3 = dYIndexP1+1;
                dDistX = (np.mod(dShadPoint.real, dShad))/dShad;
                dDistY = (np.mod(dShadPoint.imag, dShad))/dShad;
                
                dStdNormFactor = np.sqrt((1 - 2 * dDistY + 2 * (dDistY**2)) * (1 - 2 * dDistX + 2 * (dDistX**2)));
                
                dSample1 = mtShadowingSamples[int(dYIndexP1),int(dXIndexP1),7];
                dSample2 = mtShadowingSamples[int(dYIndexP2),int(dXIndexP2),7]; 
                dSample3 = mtShadowingSamples[int(dYIndexP3),int(dXIndexP3),7]; 
                dSample4 = mtShadowingSamples[int(dYIndexP4),int(dXIndexP4),7];
                
                dShadowingC = ((1-dDistY)*(dSample1*(1-dDistX)+dSample2*(dDistX)) + (dDistY*(dSample3*(1-dDistX) + dSample4*(dDistX))))/dStdNormFactor;
                
                for iMap in range(0,7):
                    dSample1 = mtShadowingSamples[int(dYIndexP1),int(dXIndexP1),iMap];
                    dSample2 = mtShadowingSamples[int(dYIndexP2),int(dXIndexP2),iMap]; 
                    dSample3 = mtShadowingSamples[int(dYIndexP3),int(dXIndexP3),iMap]; 
                    dSample4 = mtShadowingSamples[int(dYIndexP4),int(dXIndexP4),iMap];
                    
                    dShadowingERB = ((1-dDistY)*(dSample1*(1-dDistX)+dSample2*(dDistX)) + (dDistY*(dSample3*(1-dDistX) + dSample4*(dDistX))))/dStdNormFactor;
                    
                    mtShadowingCorr[il,ic,iMap] = np.sqrt(dAlphaCorr)*dShadowingC + np.sqrt(1-dAlphaCorr)*dShadowingERB;
    return mtShadowingCorr;
                    
