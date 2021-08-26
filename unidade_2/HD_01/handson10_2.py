#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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


# Parâmetros
n_bits = 100    # Número de bits
T = 50          # Tempo de símbolo
Ts = 2          # Tempo de símbolo em portadora única
K = T/Ts        # Número de subportadoras independentes
N = 2*K         # Número de pontos da IDFT

# Gerando bits aleatórios
dataIn = np.random.rand(1,n_bits)
#dataIn = np.ones((1,100))
dataIn = np.sign(dataIn-.5)

# Conversor serial paralelo
dataInMatrix = np.reshape(dataIn, (int(n_bits/4),4))

# Gerar constelação 16-QAM
seq16qam = 2*dataInMatrix[:,0]+dataInMatrix[:,1]+1j*(2*dataInMatrix[:,2]+dataInMatrix[:,3])
seq16 = seq16qam[None,:]

# Garantir propriedade da simetria
X = seq16
X = np.append(X,np.fliplr(np.conj(seq16)), axis=1)

# Construindo xn
xn = np.zeros((1,int(N)), dtype=complex)
for n in range(0,int(N)):
    for k in range(0,int(N)):
        xn[0,n] = xn[0,n] + 1/np.sqrt(N)*X[0,k]*np.exp(1j*2*np.pi*n*k/N);
        
# Construindo xt
xt = np.zeros((1, T+1), dtype=complex)
for t in range(0,T+1):
    for q in range(0,int(N)):
       xt[0,t]=xt[0,t]+1/np.sqrt(N)*X[0,q]*np.exp(1j*2*np.pi*q*t/T); 

# Plots
plt.plot(np.transpose(np.absolute(xt)),'b')
plt.stem(np.transpose(np.absolute(xn)), linefmt="r", markerfmt="ro",basefmt=" ", use_line_collection="True")
plt.title("Sinais OFDM")
plt.xlabel("Tempo")
plt.legend(["x(t)", "Xn"])
plt.show()