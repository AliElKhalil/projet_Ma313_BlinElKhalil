# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 20:12:50 2021

@author: aliel
"""

import numpy as np
import matplotlib.pyplot as plt
from Fonctions import *
import random as rdm
from math import *
from Projet_Partie_1 import *


def RegressionLineaireEllipse(X,Y):
    X=[X]
    Y=[Y]
    Xt=np.transpose(X)
    Yt=np.transpose(Y)
    n=np.shape(Xt)
    A=np.concatenate((Xt*Xt,Xt*Yt,Yt*Yt,Xt,Yt,np.ones(n)),axis=1)
    B=np.ones(n)
    M=ResolMCSVD(A,B)
    return M

 
def TraceEllipse(a,b,alpha,beta,teta):
    c=max(a,b)
    gamma=max(abs(alpha),abs(beta))
    x = -np.linspace(-c-gamma,c+gamma,1000)
    y = np.linspace(-c-gamma,c+gamma,1000)
    X,Y = np.meshgrid(x,y)
    A=(np.cos(teta)**2)/(a**2)+(np.sin(teta)**2)/(b**2)
    B=2*np.cos(teta)*np.sin(teta)
    C=(np.cos(teta)**2)/(b**2)+(np.sin(teta)**2)/(a**2)
    h=A+C
    D=-2*alpha*A-B*beta*h
    E=-2*beta*C-B*alpha*h
    F=(((alpha**2)*(np.cos(teta)**2)+(beta**2)*(np.sin(teta)**2))/(a**2))+(((beta**2)*(np.cos(teta)**2)+(alpha**2)*(np.sin(teta)**2))/(b**2))
    eqn = A*(X**2)+B*X*Y+C*(Y**2)+D*X+E*Y+F
    Z = 1
    plt.xlim([-c-gamma,c+gamma])
    plt.ylim([-c-gamma,c+gamma])
    plt.contour(X,Y,eqn,[Z])
    plt.grid()
    plt.show()
    
    
   
def TracerEllipseRegression(X,Y):
    M=RegressionLineaireEllipse(X,Y)
    A=M[0][0]
    B=M[1][0]
    C=M[2][0]
    D=M[3][0]
    E=M[4][0]
    F=M[5][0]
    teta=np.arcsin(B)/2
    b=np.sqrt((-(np.tan(teta)**2)*((np.sin(teta)**2)+(np.cos(teta)**2)))/C*(1-(A*np.tan(teta)**2)/C))
    a=np.sqrt(((b**2)*(np.cos(teta)**2))/((A*(b**2))-(np.sin(teta)**2)))
    h=A+C
    beta=(-2*E*A+B*D)/(2*A*(2*C+((h**2)*B**2)/(-2*A)))
    alpha=(D+B*beta*h)/(-2*A)
    p = -np.linspace(min(X)-1,max(X)+1,1000)
    q = np.linspace(min(Y)-1,max(Y)+1,1000)
    x,y = np.meshgrid(p,q)
    eqn = A*(x**2)+B*y*x+C*(y**2)+D*x+E*y+F
    Z = 1
    plt.scatter(X,Y)
    plt.scatter(alpha,beta)
    plt.contour(x,y,eqn,[Z])
    plt.show()
    



def RegressionLineairePolynomiale(X,Y,p):
    j=np.shape(X)[0]
    if np.shape(X)==(j,):
        X=[X]
    j=np.shape(Y)[0]
    if np.shape(Y)==(j,):
        Y=[Y]
    Xt=np.transpose(X)
    Yt=np.transpose(Y)
    M=np.ones(np.shape(Xt))
    F=np.ones(np.shape(Xt))
    for i in range(1,p+1):
        F=F*Xt
        M=np.concatenate((M,F),axis=1)
    return ResolMCSVD(M,Yt)

def ComparaisonPolynomeRegression(X,Y,p):
    x=np.linspace(min(X),max(X),10000)
    y=[]
    R=RegressionLineairePolynomiale(X, Y, p)
    for i in x:
        P=R[0][0]
        for j in range (1,p+1):
            P=P+R[j][0]*i**j
        y.append(P)
    plt.plot(x,y)
    plt.scatter(X,Y)
    
    
def ResolMCSVD(A,b):
    """
    Fonction qui résout le problème des moindres carrés en utilisant
    la décomposition en valeur singulière.

    Parameters
    ----------
    A : array
        Matrice A du système Ax=b.
    b : array
        Matrice b du système Ax=b.

    Returns
    -------
    x : array
        Solution au moindre carré par la déomposition en valeur
        singulière.

    """
    K=np.linalg.svd(A,False)
    S_p_i=np.diag(1/K[1])
    U_e=np.transpose(K[0])
    V=np.transpose(K[2])
    A_p_i=V@S_p_i@U_e
    x=np.dot(A_p_i,b)
    return x


#def RegressionLineaireEllipseV2(X,Y):
    #j=np.shape(X)[0]
    #if np.shape(X)==(j,):
    #    X=[X]
    #j=np.shape(Y)[0]
    #if np.shape(Y)==(j,):
    #    Y=[Y]
    #Xt=np.transpose(X)
    #Yt=np.transpose(Y)
    #D=np.concatenate(Xt,Yt,np.ones(np.shape(Xt)),Xt*Xt,np.sqrt(2)*Xt*Yt,Yt*Yt)
    #R=DecompositionGSGenerale(D)[1]
    #p=np.shape(R)[1]
    #R22=R[p-3:][:,p-3:]
    #w=ResolMCSVD(R22,np.zeros((1,3)))    
    #return w

def Donnees_parabole():
    X=[0.1,0.2,0.3,0.4,0.5]
    Y=[0.225,0.209,0.200,0.221,0.259]
    return X,Y