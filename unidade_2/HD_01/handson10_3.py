#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 18:39:09 2021

@author: lucasismael
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import cmath
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from itertools import cycle
import scipy.stats as st
import os
import argparse
import yaml
from random import randint   
import itertools  

# Parâmetros
n_bits = 1000    # Número de bits
T = 500          # Tempo de símbolo
Ts = 2          # Tempo de símbolo em portadora única
K = int(T/Ts)        # Número de subportadoras independentes
N = int(2*K)         # Número de pontos da IDFT
sigmas = [0, 0.1, 1] #Vetor de variâncias do ruído
#sigmas = [0]

# Gerar bits aleatórios
dataIn = np.random.rand(1,n_bits)
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
        
# Loop de variâncias
for ik in sigmas:
    # Adição de ruído
    variance = ik
    noise = np.sqrt(variance)*np.random.randn(1,N)+1j*np.sqrt(variance)*np.random.randn(1,N)
    
    # Sinal recebido = xn + ruído
    rn = xn + noise
    # DFT de rn
    Y = np.zeros((1,K), dtype="complex")
    for k in range(0,K):
        for n in range(0,N):
            Y[0, k] = Y[0, k] + 1/np.sqrt(N)*rn[0,n]*np.exp(-1j*2*np.pi*k*n/N)
    
    plt.figure()
    plt.scatter(np.real(Y), np.imag(Y), c='b', marker='.')
    plt.scatter(np.real(seq16), np.imag(seq16), c='r', marker='x')
    plt.title("Sinal com ruído de variância " + str(variance))
    plt.xlabel("In-Phase")
    plt.ylabel("Quadrature")
    plt.show()
    
    # Demodulação
    Z = np.zeros((1,K), dtype="complex")
    for k in range(0,Y.shape[1]):
        if np.real(Y[0,k]) > 0:
            if np.real(Y[0,k]) > 2:
                Z[0,k] = 3
            else:
                Z[0,k] = 1
        else: 
            if np.real(Y[0,k]) < -2:
                Z[0,k] = -3
            else:
                Z[0,k] = -1
        if np.imag(Y[0,k]) > 0:
            if np.imag(Y[0,k]) > 2:
                Z[0,k] = Z[0,k] + 1j*3
            else:
                Z[0,k] = Z[0,k] + 1j
        else:
            if np.imag(Y[0,k]) < -2:
                Z[0,k] = Z[0,k] - 1j*3
            else:
                Z[0,k] = Z[0,k] - 1j
    
    #error = len(np.nonzero(Z - X[0,:250]))
    error = np.count_nonzero(Z - X[0,:250])
    print("Para a variância de " + str(variance) + " houve " + str(error) + " símbolos errados." )
        