# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 20:57:28 2020

@author: aliel
"""

import numpy as np
import matplotlib.pyplot as plt
from Fonctions import *
import random as rdm
from math import *
from Projet_Partie_1 import *


def RegressionLineaireCercle(X,Y):
    """

    Parameters
    ----------
    X : List
        Coordonnees en abscisse des points a etudier.
    Y : List
        Coordonnees en ordonee des points a etudier.

    Returns
    -------
    M : list
        liste contenant la coordonnee en abscisse et en ordonnee du
        centre et rayon du cercle de regression lineaire.

    """
    j=np.shape(X)[0]
    if np.shape(X)==(j,):
        X=[X]
    j=np.shape(Y)[0]
    if np.shape(Y)==(j,):
        Y=[Y]
    Xt=np.transpose(X)
    Yt=np.transpose(Y)
    n=np.shape(Xt)
    A=2*np.concatenate((Xt,Yt,1/2*np.ones(n)),axis=1)
    B=Xt*Xt+Yt*Yt
    M=ResolMCQR(A,B)
    return M


def CourbePoint(X,Y):
    """
    Fonction qui trace le nuage de point à partir des coordonnées

    Parameters
    ----------
    X : List
        Coordonnees en abscisse des points a etudier.
    Y : List
        Coordonnees en ordonee des points a etudier.

    Returns
    -------
    None.

    """
    plt.scatter(X,Y)
    plt.show()

def TraceCercle(a,b,r):
    """
    Fonction qui trace un cercle, avec un centre et un rayon donnés.

    Parameters
    ----------
    a : float
        Coordonnée en abscisse du centre du cercle.
    b : float
        Coordonnée en ordonnée du centre du cercle.
    r : float
        Rayon du cercle.

    Returns
    -------
    None.

    """
    theta=np.linspace(0,2*np.pi,1000)
    x=r*np.cos(theta)+a
    y=r*np.sin(theta)+b
    plt.plot(x,y)
    plt.show()

def ComparaisonPointsRegression(X,Y):
    """
    Fonction qui trace le nuage de points et la regressions linéaire
    d'une liste de points

    Parameters
    ----------
    X : List
        Coordonnees en abscisse des points a etudier.
    Y : List
        Coordonnees en ordonnee des points a etudier.

    Returns
    -------
    None.

    """
    [[a],[b],[g]]=RegressionLineaireCercle(X,Y)
    R=np.sqrt(g+a**2+b**2)
    alpha=a
    beta=b
    theta=np.linspace(0,2*np.pi,1000)
    x=R*np.cos(theta)+a
    y=R*np.sin(theta)+b
    plt.plot(x,y,label="régression linéaire")
    plt.scatter(a,b,label="centre du cercle")
    plt.scatter(X,Y,label="points étudiés")
    plt.title("Comparaison des points avec la regression linéaire")
    plt.legend()
    plt.legend(loc='upper left')
    plt.show()
    

def ResolutionPartie3():
    X,Y=donnees_partie3()
    ComparaisonPointsRegression(X, Y)
    
    
def RegressionLineaireEllipse(X,Y):
    X=[X]
    Y=[Y]
    Xt=np.transpose(X)
    Yt=np.transpose(Y)
    n=np.shape(Xt)
    A=np.concatenate((Xt*Xt,Xt*Yt,Yt*Yt,Xt,Yt,np.ones(n)),axis=1)
    B=np.ones(n)
    M=ResolMCNP(A,B)
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
    plt.ylim([-c-gamma*c,c+gamma])
    plt.contour(X,Y,eqn,[Z])
    plt.grid()
    plt.show()
    
    
   
def tracerEllipseRegression(X,Y):
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
    beta=b*np.sqrt(E)
    if a-b>0:
        c=a
    else:
        c=b
    if np.abs(alpha)-np.abs(beta)>0:
        gamma=np.abs(alpha)
    else:
        gamma=np.abs(beta)
    p = -np.linspace(-c-gamma,c+gamma,1000)
    q = np.linspace(-c-gamma,c+gamma,1000)
    x,y = np.meshgrid(p,q)
    eqn = A*(x**2)+B*y*x+C*(y**2)+D*x+E*y+F
    Z = 1
    plt.contour(x,y,eqn,[Z])
    plt.show()
    



def RegressionLineairePolynomial(X,Y,p):
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
        M=np.concatenate(F,M,axis=1)
    return ResolMCQR(M,Yt)

def TracerPolynomeRegression(X,Y,p):
    x=np.linspace(min(X),max(X),10000)
    y=[]
    R=RegressionLineairePolynomial(X, Y, p)
    for i in x:
        P=0
        for j in range (0,p+1):
            P=P+R[j]**j
        y.append(P)
    plt.plot(x,y)
    plt.plot(X,Y)
    
    
def ResolMCSVD(A,b):
    S_p_i=1/np.linalg.svd(A)[1]
    U_e=np.matrix.H(np.linalg.svd(A)[0])
    V=np.matrix.H(np.linalg.svd(A)[2])
    A_p_i=np.dot(V,S_p_i,U_e)
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



    
            
    

    