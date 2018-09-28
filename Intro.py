# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 18:44:11 2018

@author: Admin
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize
#%%
def ldata(archive):
    f=open(archive)
    data=[]
    for line in f:
        line=line.strip()
        col=line.split()
        data.append(col)
    return data
#%%
prot = ldata('C:/Users/Admin/Documents/GitHub/Redes_TP2/data/yeast_AP-MS.txt')
bina = ldata('C:/Users/Admin/Documents/GitHub/Redes_TP2/data/yeast_Y2H.txt')
lit = ldata('C:/Users/Admin/Documents/GitHub/Redes_TP2/data/yeast_LIT.txt')
litreg = ldata('C:/Users/Admin/Documents/GitHub/Redes_TP2/data/yeast_LIT_Reguly.txt')
essential = ldata('C:/Users/Admin/Documents/GitHub/Redes_TP2/data/Essential_ORFs_paperHe.txt')
#%% LIMPIO LOS DATOS
ess = [a[1] for a in essential[2:1158]]
litr = [a[0:2] for a in litreg[1:]]
#%%
P = nx.Graph()
P.add_edges_from(prot)

B = nx.Graph()
B.add_edges_from(bina)

L = nx.Graph()
L.add_edges_from(lit)

LR = nx.Graph()
LR.add_edges_from(litreg)
#%%
def topologia(redes): #dar vector con redes ordenadas
    G = redes
    N = np.empty_like(G)          # número de nodos de la red
    E = np.empty_like(G)          # número de enlaces de la red
    G_mean = np.empty_like(G)     # grado medio de la red
    G_max = np.empty_like(G)      # grado máximo
    G_min  = np.empty_like(G)     # grado mínimo
    dens = np.empty_like(G)       # densidad de la red
    clu = np.empty_like(G)        # coeficiente de clustering local
    clu_d = np.empty_like(G)      # coeficiente de clustering global/transitividad
    diam = np.empty_like(G)       # diámetro de la red
    
    def directed(A):
    for i in range(0,len(A)):
        for r in range(i+1,len(A)):
            if A[i][0] == A[r][1] and A[i][1] == A[r][0]:
                return "SI"
    return "NO"