#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 18:26:05 2018

@author: camilo
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
#
def ldata(archive):
    f=open(archive)
    data=[]
    for line in f:
        line=line.strip()
        col=line.split()
        data.append(col)
    return data

APMS_data = ldata('yeast_AP-MS.txt')
LIT_data = ldata('yeast_LIT.txt')
LIT_Reguly_data = ldata('yeast_LIT_Reguly.txt')
Y2H_data = ldata('yeast_Y2H.txt')
ess_data = ldata('Essential_ORFs_paperHe.txt')

#Primero me interesa ver las propiedaddes estructurales de las redes de 
#interacción de proteínas

#Genero los grafost
yeastY2H = nx.Graph(Y2H_data)
yeastLIT = nx.Graph(LIT_data)
yeastAPMS = nx.Graph(APMS_data)
yeastLIT_Reguly = nx.Graph(LIT_Reguly_data)

#Tengo que hacer reguly_data como una lista de las interacciones entre proteínas
for h in range(1,21135):
    Reguly_data = [list(),list()]
    Reguly_data[h] = [LIT_Reguly_data[h][0],LIT_Reguly_data[h][1]]


#Genero una función que me devuelva el overlap de edges entre dos grafos
def overlap(A,B):
    m = 0
    if len(A)>len(B):
        l_aux = len(A)-len(B)
        B1 = list()
        for r in range(0,l_aux):
            B1[r] = ['hello world','hello world']#Algoo para subsanar
        #Más facil es ver cuantas duplas se repiten por lista, compara A con B
        B2  = B.append(B1)
        A2 = A
    if len(B)>len(A):
        l_aux = len(B)-len(A)
        A1 = list()
        for l in range(0,l_aux):
            A1[l] = ['steve','madden']
            
        A2 = A.append(A1)
        B2 = B
    for i in range(0,len(A2)):
        for j in range(0,len(A2)):
            if A2[i][0] == B2[j][0] and A2[i][1] == B2[j][1]:
                m = m + 1
            if A2[i][0] == B2[j][1] and A2[i][1] == B2[j][0]:
                m = m + 1
    return m/len(A2)

def overlap1(A,B):
    lst3 = [value for value in A if value in B]
    return len(lst3)/float(len(A))

#def overlap2(A,B):
#    m = 0
#    for i in range(0,len(A)):
#        for j in range(0,len(B)):
#            if A[i][0] == B[j][0] and A[i][1] == B[j][1]:
#                m = m + 1
#            if A[i][0] == B[j][1] and A[i][1] == B[j][0]:
#                m = m + 1
#    return m/len(A)
    
def overlap2(A,B):
    m = 0
    for nodo_a in A:
        for nodo_b in B:
            flag = 0
            if nodo_a[0] == nodo_b[0] and nodo_a[1] == nodo_b[1]:
                flag = 1
            if nodo_a[0] == nodo_b[1] and nodo_a[1] == nodo_b[0]:
                flag = 1
            m = m + flag
    return m/float(len(A))

def algo(A):
    for enlaces in A:
        if len(enlaces) != 2:
            print('cuidado')
            break
            

#for nodo_a in A:
#    
#    if nodo_a[0]

ess = [a[1] for a in ess_data[2:1158]]
litr = [a[0:2] for a in LIT_Reguly_data[1:]]

def lista(A):
    C = []
    for nodo_a in A:
        C.append(nodo_a[0]+' '+nodo_a[1])
    return C    

test_reg = lista(litr)
test_apms = lista(APMS_data)
test_y2h = lista(Y2H_data)
test_lit = lista(LIT_data)

def arreglo(A):
    C = []
    for nodo_a in A:
        C.append(nodo_a[0]+' '+nodo_a[1])    
    arreglado = set(C)
    D=[]
    for i in rqange(0, len(list(arreglado))):
        D[i] = list(arreglado)[i].split()
    for r in range(0,len(D)):
        for k in range(0,len(D)):
            if D[r][0]==D[k][1] and D[r][1] == D[k][0]:
                D.pop(k)
            
    return D

        
#arreglo limpmia las tuplas de las listas que estan repetidas con set, y despues queremos borrar tambien las tuplas que estan dadaas vueltq.  Todo esto para contar bien el overlap (solo interesa una vez cada interaccion)
    
ovlp_Y2H_w_LIT = overlap2(Y2H_data,LIT_data) #fracción de interacciones en Y2H que aparecen en LIT
ovlp_Y2H_w_APMS = overlap2(Y2H_data,APMS_data)
ovlp_Y2H_w_reguly = overlap2(Y2H_data,LIT_Reguly_data)
#ovlp_Y2H_w_ess = overlap2(Y2H_data,ess_data)

ovlp_LIT_w_Y2H = overlap2(LIT_data,Y2H_data) #fracción de interacciones en LIT que aparecen en Y2H
ovlp_LIT_w_APMS = overlap2(LIT_data,APMS_data) #f
ovlp_LIT_w_reguly = overlap2(LIT_data,LIT_Reguly_data)
#ovlp_LIT_w_ess = overlap2(LIT_data,ess_data)

ovlp_APMS_w_Y2H = overlap2(APMS_data,Y2H_data)
ovlp_APMS_w_LIT = overlap2(APMS_data,LIT_data)
ovlp_APMS_w_reguly = overlap2(APMS_data,LIT_Reguly_data)
#ovlp_APMS_w_ess = overlap2(APMS_data,Y2H_data)

ovlp_reguly_w_Y2H = overlap2(LIT_Reguly_data,Y2H_data)
ovlp_reguly_w_LIT = overlap2(LIT_Reguly_data,LIT_data)
ovlp_reguly_w_APMS = overlap2(LIT_Reguly_data,APMS_data)
#ovlp_reguly_w_ess = overlap2(LIT_Reguly_data, ess_data)

ovlp_ess_w_Y2H = overlap2(LIT_Reguly_data,Y2H_data)
ovlp_ess_w_LIT = overlap2(LIT_Reguly_data,LIT_data)
ovlp_ess_w_APMS = overlap2(LIT_Reguly_data,APMS_data)
#ovlp_ess_w_reguly = overlap2(LIT_Reguly_data, ess_data)

firstcol = ['Y2H',ovlp_Y2H_w_LIT,ovlp_Y2H_w_APMS,ovlp_Y2H_w_reguly];
secondcol = [ovlp_LIT_w_APMS,'LIT',ovlp_LIT_w_APMS, ovlp_LIT_w_reguly];
thirdcol = [ovlp_APMS_w_Y2H, ovlp_APMS_w_LIT, 'APMS', ovlp_APMS_w_reguly];
fourthcol = [ovlp_reguly_w_Y2H,ovlp_reguly_w_LIT,ovlp_reguly_w_APMS,'Reguly'];
#fifthcol = [ovlp_ess_w_Y2H,ovlp_ess_w_LIT,ovlp_ess_w_APMS, ovlp_ess_w_reguly,'Essentials'];

tabla_overlaps = pd.DataFrame({"1":firstcol,"2":secondcol,"3":thirdcol,"4":fourthcol})





