"""
Created on Tue Dec  8 17:16:36 2020

@author: aliel
"""

"""
Created on Tue Dec  8 17:16:36 2020

@author: aliel
"""

import numpy as np
from Projet_Partie_1 import *
from Fonctions import *

def SystemeExercice1():
    A1=np.array([[1,2],[2,3],[-1,2]])
    b1=np.array([[12],[17],[6]])
    return A1,b1

def SystemeExercice2():
    A2=np.array([[1,21],[-1,-5],[1,17],[1,17]])
    b2=np.array([[3],[-1],[1],[1]])
    return A2,b2

def SystemeExercice3():
    x=[0.3,-2.7,-1.9,1.2,-2.6,2.7,2.0,-1.6,-0.5,-2.4]
    y=[2.8,-9.4,-4.5,3.8,-8.0,3.0,3.9,-3.5,1.3,-7.6]
    A=np.zeros((len(x),3))
    Y=np.zeros((len(x),1))
    for i in (0,len(x)):
        A[i]=[1,x[i],x[i]**2]
        Y[i]=y[i]
    return A,y


def ResultatExercice(k): 
    """
    Fonction qui donne les resultats des exercices 1, 2 et 3 
    selon les 3 méthodes definit partie 1, à savoir le systeme 
    d'equation normale, la decomposition QR et np.linalg.lstsq().
    Elle rend donc (ResolMCEN(A,b),ResolMCQR(A,b),ResolMCNP(A,b))
    """
    
    #k=input("verifier avec le système de l'exercice =?")
    if k ==1:
        A,b=SystemeExercice1()
        return ResolMCEN(A,b),ResolMCQR(A,b),ResolMCNP(A,b)
    elif k==2:
        A,b=SystemeExercice2()
        return ResolMCEN(A,b),ResolMCQR(A,b),ResolMCNP(A,b)
    elif k==3:
        A,b= SystemeExercie3()
        return ResolMCEN(A,b),ResolMCQR(A,b),ResolMCNP(A,b)
    

