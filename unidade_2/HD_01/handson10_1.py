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

# Passo 1: Geração de fases aleatórias
phi_k = 2*np.pi*np.random.rand()
phi_j = 2*np.pi*np.random.rand()

# Passo 2: Geração de sinais amostrados
M = 50;
m = np.arange(0,M)
x_k = np.sin(4*np.pi*m/5+phi_k)
n = 1
x_j_1 = np.sin(4*np.pi*m/5+2*np.pi*m*n/M+phi_j)
n = 2
x_j_2 = np.sin(4*np.pi*m/5+2*np.pi*m*n/M+phi_j)
n = 3 
x_j_3 = np.sin(4*np.pi*m/5+2*np.pi*m*n/M+phi_j)

# Passo 3: Verificação de ortogonalidade
Sum1 = np.sum(x_k*x_j_1)
print('O resultado para n=1 é: '+ str(Sum1))
Sum2 = np.sum(x_k*x_j_2)
print('O resultado para n=2 é: '+ str(Sum2))
Sum3 = np.sum(x_k*x_j_3)
print('O resultado para n=3 é: '+ str(Sum3))