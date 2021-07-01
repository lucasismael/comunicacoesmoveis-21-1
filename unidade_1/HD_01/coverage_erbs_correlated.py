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
from fCorrShadowing import fCorrShadowing 

def fDrawSector (dR, dCenter):
    vtHex=np.zeros([1,0])
    for ie in range(1,7):
        vtHex = np.append(vtHex, dR*complex(np.cos((ie-1)*np.pi/3), np.sin((ie-1)*np.pi/3)))
    vtHex = vtHex+dCenter
    vtHexp = np.append(vtHex, vtHex[0])
    #print(vtHex)
    plt.plot(vtHexp.real,vtHexp.imag, 'k')
    #plt.grid()
    #plt.tight_layout()
    #plt.show()
    
def fDrawDeploy(dR, vtBs):
    for iBsd in range(len(vtBs)):
        fDrawSector(dR, vtBs[iBsd])
    plt.plot(vtBs.real,vtBs.imag, 'sk')
    plt.axis('scaled')
    #plt.grid()
    plt.tight_layout()
    plt.show()
    

# Entrada de parâmetros
dR = 200; #Raio do Hexágono
dFc = 800;
dShad = 50;
dSigmaShad = 8;
dAlphaCorr = 0.5;

dPasso = 10;
dRmin = dPasso;
dIntersiteDistance = 2*np.sqrt(3/4)*dR;

dDimXOri = 5*dR;
dDimYOri = 6*np.sqrt(3/4)*dR;
dPtdBm = 57;
dPtLinear = 10**(dPtdBm/10)*(1e-3);
dHMob = 5;
dHBs = 30;
dAhm = 3.2*(np.log10(11.75*dHMob))**2 - 4.97;

vtBs= [0];
dOffset = np.pi/6;
for iBs in range(2,8):
    vtBs = np.append(vtBs, dR*np.sqrt(3)*np.exp(complex(0,(iBs-2)*np.pi/3 + dOffset)));
vtBs = vtBs + complex(dDimXOri/2,dDimYOri/2);

dDimY = dDimYOri + np.mod(dDimYOri, dPasso);
dDimX = dDimXOri + np.mod(dDimXOri, dPasso);
[mtPosx, mtPosy] = np.meshgrid(np.arange(0,dDimXOri+1,dPasso), np.arange(0,dDimYOri+1,dPasso));
mtPontosMedicao = mtPosx + 1j*mtPosy; 

mtPosEachBS = np.empty([mtPosx.shape[0],mtPosx.shape[1],len(vtBs)], dtype=complex);
mtPowerEachBSdBm = np.empty([mtPosx.shape[0],mtPosx.shape[1],len(vtBs)]);
mtPowerEachBSShaddBm = np.empty([mtPosx.shape[0],mtPosx.shape[1],len(vtBs)]);
mtPowerEachBSShadCorrdBm = np.empty([mtPosx.shape[0],mtPosx.shape[1],len(vtBs)]);
mtPowerFinaldBm = -np.inf*np.ones([mtPosy.shape[0],mtPosy.shape[1]]);
mtPowerFinalShaddBm = -np.inf*np.ones([mtPosy.shape[0],mtPosy.shape[1]]);
mtPowerFinalShadCorrdBm = -np.inf*np.ones([mtPosy.shape[0],mtPosy.shape[1]]);

for x in range(mtPosx.shape[0]):
    for y in range (mtPosx.shape[1]):
        mtPosEachBS[x,y,:] = complex(mtPosx[x,y],mtPosy[x,y]);
        
mtShadowingCorr = fCorrShadowing(mtPontosMedicao, dShad, dAlphaCorr, dSigmaShad, dDimXOri, dDimYOri);
meanstd = [];
for z in range(0,7):
    stdproof = np.std(mtShadowingCorr[:,:, z]);
    meanstd = np.append(meanstd, stdproof);
    print('Desvio-padrão da ERB ' + str(z) + ': ' + str(stdproof));
meanstdproof = np.mean(meanstd);
# if meanstdproof < dSigmaShad:
#     meanstdproof = np.ceil(meanstdproof);
# elif meanstdproof > dSigmaShad:
#     meanstdproof = np.floor(meanstdproof);
print('Desvio-padrão médio obtido: ' + str(meanstdproof));
                       
for iBsD in range(len(vtBs)):
    mtPosEachBS[:,:,iBsD] = mtPosEachBS[:,:,iBsD]-vtBs[iBsD];
    mtDistEachBS = np.abs(mtPosEachBS[:,:,iBsD]);
    mtDistEachBS[mtDistEachBS < dRmin] = dRmin;
    # Okumura-Hata (cidade urbana) - dB
    mtPldB = 69.55 + 26.16*np.log10(dFc) + (44.9 - 6.55*np.log10(dHBs))*np.log10(mtDistEachBS/1e3) - 13.82*np.log10(dHBs) - dAhm;
    
    mtShadowing = dSigmaShad*np.random.randn(mtPosy.shape[0], mtPosy.shape[1]);
    
    mtPowerEachBSdBm[:,:,iBsD] = dPtdBm - mtPldB;
    
    mtPowerEachBSShaddBm[:,:,iBsD] = dPtdBm - mtPldB + mtShadowing;
    
    mtPowerEachBSShadCorrdBm[:,:, iBsD] = dPtdBm - mtPldB + mtShadowingCorr[:,:,iBsD];
    #plt.figure(iBsD);
    #plt.plot(mtPosEachBS[:,:,iBsD].real,mtPosEachBS[:,:,iBsD].imag, 'bo');
    #plt.pcolormesh(mtPosx, mtPosy, mtPowerEachBSdBm[:,:,iBsD], cmap='jet');
    #plt.colorbar();
    #plt.title("ERB" + str(iBsD));
    #fDrawDeploy(dR, vtBs);
    #plt.axis('scaled')
    mtPowerFinaldBm = np.maximum(mtPowerFinaldBm, mtPowerEachBSdBm[:,:,iBsD]);
    mtPowerFinalShaddBm = np.maximum(mtPowerFinalShaddBm, mtPowerEachBSShaddBm[:,:,iBsD]);
    mtPowerFinalShadCorrdBm = np.maximum(mtPowerFinalShadCorrdBm, mtPowerEachBSShadCorrdBm[:,:,iBsD]);
    #plt.show();
plt.figure();
plt.pcolormesh(mtPosx,mtPosy,mtPowerFinaldBm, cmap='hsv');
plt.colorbar();
fDrawDeploy(dR, vtBs);
plt.axis('scaled');
#plt.tight_layout()
plt.title("REM MAP - Sem Shadowing");
plt.figure();
plt.pcolormesh(mtPosx,mtPosy,mtPowerFinalShaddBm, cmap='hsv');
plt.colorbar();
fDrawDeploy(dR, vtBs);
plt.axis('scaled');
#plt.tight_layout()
plt.title("REM MAP - Com Shadowing");
plt.figure();
plt.pcolormesh(mtPosx,mtPosy,mtPowerFinalShadCorrdBm, cmap='hsv');
plt.colorbar();
fDrawDeploy(dR, vtBs);
plt.axis('scaled');
#lt.tight_layout()
plt.title("REM MAP - Com Shadowing Correlacionado");
plt.show();




    
