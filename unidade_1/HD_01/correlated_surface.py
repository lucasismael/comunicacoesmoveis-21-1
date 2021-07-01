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
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from itertools import cycle
import scipy.stats as st
import os
import argparse
import yaml
import matplotlib
from random import randint   
import itertools  

dR = 150;
dShad = 50;
dPasso = 7;
dSigmaShad = 8;

dDimXOri = 5*dR;
dDimYOri = 6*np.sqrt(3/4)*dR;

dDimY = np.ceil(dDimYOri + np.mod(dDimYOri, dPasso));
dDimX = np.ceil(dDimXOri + np.mod(dDimXOri, dPasso));
[mtPosx, mtPosy] = np.meshgrid(np.arange(0,dDimX+1,dPasso), np.arange(0,dDimY+1,dPasso));
mtPontosMedicao = mtPosx + 1j*mtPosy;


dDimYS = np.ceil(dDimYOri + np.mod(dDimYOri, dShad));
dDimXS = np.ceil(dDimXOri + np.mod(dDimXOri, dShad));
[mtPosxshad, mtPosyshad] = np.meshgrid(np.arange(0,dDimXS+1,dShad), np.arange(0,dDimYS+1,dShad));
mtPosShad = mtPosxshad + 1j*mtPosyshad;

mtShadowingSamples = dSigmaShad*np.random.randn(mtPosy.shape[0], mtPosy.shape[1]);
mtShadowingCorr = np.empty([mtPosy.shape[0], mtPosy.shape[1]]);

dSizel = mtPontosMedicao.shape[0];
dSizec = mtPontosMedicao.shape[1];

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

for il in range(0,dSizel):
    for ic in range(0,dSizec):
        
        dshadPoint = mtPontosMedicao[il,ic];
        
        dXIndexP1 = dshadPoint.real/dShad;
        dYIndexP1 = dshadPoint.imag/dShad;
        
        if (np.mod(dXIndexP1,1) == 0 and np.mod(dYIndexP1,1) == 0):
            dXIndexP1 = np.floor(dXIndexP1);
            dYIndexP1 = np.floor(dYIndexP1);
            # plt.plot(mtPosShad.real,mtPosShad.imag, 'ko');
            # plt.plot(complex(mtPosShad[int(dYIndexP1),int(dXIndexP1)]), 'g*');
            # plt.axis('equal');
            # plt.xlim([-2*dShad+mtPosShad[int(dYIndexP1),int(dXIndexP1)].real,2*dShad+mtPosShad[int(dYIndexP1),int(dXIndexP1)].real])
            # plt.ylim([-2*dShad+mtPosShad[int(dYIndexP1),int(dXIndexP1)].imag,2*dShad+mtPosShad[int(dYIndexP1),int(dXIndexP1)].imag]);
            # plt.show();
            # print('O ponto de medição é um ponto de grade');
            mtShadowingCorr[il,ic] = mtShadowingSamples[int(dYIndexP1),int(dXIndexP1)];
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
            # plt.plot(mtPosShad.real,mtPosShad.imag, 'ko');
            # plt.plot(dshadPoint.real, dshadPoint.imag, 'sr');
            # plt.plot(mtPosShad[int(dYIndexP1),int(dXIndexP1)].real, mtPosShad[int(dYIndexP1),int(dXIndexP1)].imag, 'b*');
            # plt.plot(mtPosShad[int(dYIndexP2),int(dXIndexP2)].real, mtPosShad[int(dYIndexP2),int(dXIndexP2)].imag, 'b*');
            # plt.plot(mtPosShad[int(dYIndexP3),int(dXIndexP3)].real, mtPosShad[int(dYIndexP3),int(dXIndexP3)].imag, 'b*');
            # plt.plot(mtPosShad[int(dYIndexP4),int(dXIndexP4)].real, mtPosShad[int(dYIndexP4),int(dXIndexP4)].imag, 'b*');
            # plt.axis('equal');
            # plt.xlim([-2*dShad+mtPosShad[int(dYIndexP3),int(dXIndexP3)].real,2*dShad+mtPosShad[int(dYIndexP4),int(dXIndexP4)].real])
            # plt.ylim([-2*dShad+mtPosShad[int(dYIndexP3),int(dXIndexP3)].imag,2*dShad+mtPosShad[int(dYIndexP1),int(dXIndexP1)].imag]);
            # plt.show();
        
            dDistX = (np.mod(dshadPoint.real, dShad))/dShad;
            dDistY = (np.mod(dshadPoint.imag, dShad))/dShad;
            #print('X = ' + str(dDistX) + ' e Y = ' + str(dDistY));
            
            dStdNormFactor = np.sqrt((1 - 2 * dDistY + 2 * (dDistY**2)) * (1 - 2 * dDistX + 2 * (dDistX**2)));
            
            dSample1 = mtShadowingSamples[int(dYIndexP1),int(dXIndexP1)];
            dSample2 = mtShadowingSamples[int(dYIndexP2),int(dXIndexP2)]; 
            dSample3 = mtShadowingSamples[int(dYIndexP3),int(dXIndexP3)]; 
            dSample4 = mtShadowingSamples[int(dYIndexP4),int(dXIndexP4)]; 
            
            mtShadowingCorr[il,ic] = ((1-dDistY)*(dSample1*(1-dDistX)+dSample2*(dDistX)) + (dDistY*(dSample3*(1-dDistX) + dSample4*(dDistX))))/dStdNormFactor;
            
surf = ax.plot_surface(mtPosx, mtPosy, mtShadowingCorr, cmap = 'jet');
fig.colorbar(surf, shrink=0.5, aspect=5);
plt.show();