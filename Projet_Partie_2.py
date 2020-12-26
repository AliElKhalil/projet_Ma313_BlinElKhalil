"""
Created on Tue Dec  8 17:16:36 2020

@author: aliel
"""

import numpy as np
from Projet_Partie_1 import *
from Fonctions import *
import random as rdm
from math import *

def SystemeExercice(k):
    if k==1:
        A1=np.array([[1,2],[2,3],[-1,2]])
        b1=np.array([[12],[17],[6]])
        return A1,b1
    if k==2:
        A2=np.array([[1,21],[-1,-5],[1,17],[1,17]])
        b2=np.array([[3],[-1],[1],[1]])
        return A2,b2
    if k==3:
        x=[0.3,-2.7,-1.9,1.2,-2.6,2.7,2.0,-1.6,-0.5,-2.4]
        y=[2.8,-9.4,-4.5,3.8,-8.0,3.0,3.9,-3.5,1.3,-7.6]
        A=np.zeros((len(x),3))
        Y=np.zeros((len(x),1))
        for i in range (0,len(x)):
            A[i]=[1,x[i],x[i]**2]
            Y[i]=y[i]
        return A,y


def ResultatExercice(k):
    """

    Parameters
    ----------
    k : Integer
        Numéro de l'exercice duquel on veut le résultat du système.

    Returns
    -------
    tuple
        résultat du problème des moindres carrés selon la méthode d'équation
        normale, avec une décomposition QR et avec la fonction linalg.lstsq
        de numpy.

    """
    A,b=SystemeExercice(k)
    return ResolMCEN(A,b),ResolMCQR(A,b),ResolMCNP(A,b)



def VerificationMinimum():
    G=[]
    for k in [1,2,3]:
        for i in [0,1,2]:
            c=True
            X=ResultatExercice(k)[i]
            A,b=SystemeExercice(k)
            N=np.linalg.norm(np.dot(A,X)-b)
            p=np.shape(X)[0]
            e=(10**-3)/(np.sqrt(p))
            j=0
            while j<=10**6 and c==True:
                x=X+rdm.uniform(0,e)
                j+=1
                n=np.linalg.norm(np.dot(A,x)-b)
                if n-N<0:
                    c=False
                    G.append((k,i,x,n-N))
    if G==[]:
        return True
    else :
        return False,G
