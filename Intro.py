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
LR.add_edges_from(litr)
#%%
def directed(A):
    for i in range(0,len(A)):
        for r in range(i+1,len(A)):
            if A[i][0] == A[r][1] and A[i][1] == A[r][0]:
                return "SI"
    return "NO"

def topologia(redes, edges): #dar vector con redes ordenadas y sus bases de datos
    
    R = redes
    L = edges
    
    N = np.empty_like(R)          # número de nodos de la red
    E = np.empty_like(R)          # número de enlaces de la red
    R_mean = np.empty_like(R)     # grado medio de la red
    R_max = np.empty_like(R)      # grado máximo
    R_min  = np.empty_like(R)     # grado mínimo
    dens = np.empty_like(R)       # densidad de la red
    clu = np.empty_like(R)        # coeficiente de clustering local
    clu_d = np.empty_like(R)      # coeficiente de clustering global/transitividad
    diam = np.empty_like(R)       # diámetro de la red
    
    for i in range(len(R)):
        N[i] = R[i].number_of_nodes()
        E[i] = R[i].number_of_edges()
        R_mean[i] = np.mean([a[1] for a in list(R[i].degree())])
        R_max[i] = np.max([a[1] for a in list(R[i].degree())])
        R_min[i] = np.min([a[1] for a in list(R[i].degree())])
        dens[i] = nx.density(R[i])
        clu[i] = nx.average_clustering(R[i])
        clu_d[i] = nx.transitivity(R[i])
        diam[i] = nx.diameter(max(nx.connected_component_subgraphs(R[i]), key=len))
    
    dirigido = np.empty_like(R)
    for i in range(len(dirigido)):
        dirigido[i] = directed(L[i])
        
    return N, E, R_mean, R_max, R_min, dens, clu, clu_d, diam, dirigido
#%%
A = topologia([P,B,L,LR], [prot, bina, lit, litr])
#%%
tabla = pd.DataFrame({"Red":["Binarias","Proteicas","Literatura", "Literatura Regulada"],"# de nodos":A[0],"# total de enlaces":A[1],"Grado medio":A[2],"Coef. de Clust. red": A[6],"Dirigida?":A[9]})
print(tabla)
#%% #función para asignar si un nodo es esenecial o no para todas las redes al mismo tiempo
def daresencialidad(redes, ess):
    R = redes
    esencial = []#np.empty_like((len(R), 1)) ##vamos a tener que usar un diccionario porque hay distinto numero de nodos
    for i in range(len(R)):
        D = dict()
        a = list(R[i].nodes())
        for j in range(len(a)):
            D[j] = 0
            for l in range(len(ess)):
                if a[j] == ess[l]:
                    D[j] = 1
                    break
        esencial.append(D)
    
    return esencial
#%% #meto el atributo esencialidad en los nodos de cada red
def atribuir(redes, ess):
    R = redes
    for i in range(len(R)):
        for n,e,g in zip(R[i].nodes, daresencialidad(redes, ess)[i].values(), [a[1] for a in list(R[i].degree())]):
            R[i].nodes[n]['Esencialidad'] = e
            R[i].nodes[n]['Grado'] = g
            
#%%
atribuir([P,B,L,LR],ess)
#%%
def eshub(redesconatributos, porcentajes): #np.linspace(0,1,11)
    ess_percent = []  #aca voy a appendear los vectores para cada red, osea que va a tener len = len(R)
    R = redesconatributos
    for i in range(len(R)):
        a = []
        for j,val in enumerate(porcentajes):
            cant_nodos_esenciales = np.sum([a[1]['Esencialidad'] for a in list(sorted(R[i].nodes.data(),key = lambda x: -x[1]['Grado']))[0:int(val*R[i].number_of_nodes())]])
            a.append(cant_nodos_esenciales/int(val*R[i].number_of_nodes())) #guardo el porcentaje de nodos esenciales hasta el número de nodos considerados como hubs por vez
        ess_percent.append(a)
    return ess_percent

#%%
#ess percent va a ser un vector con tantos vectores como redes
plt.figure(1)

plt.plot(porcentajes, ess_percent[0], label = '%s'%redesconatributos[0])
plt.plot(porcentajes, ess_percent[1], label = '%s'%redesconatributos[1])   
plt.plot(porcentajes, ess_percent[2], label = '%s'%redesconatributos[2])   
plt.plot(porcentajes, ess_percent[3], label = '%s'%redesconatributos[3])   
plt.grid(True)
plt.legend()             
