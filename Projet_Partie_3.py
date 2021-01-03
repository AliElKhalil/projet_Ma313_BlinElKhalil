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
from Projet_Partie_4 import *

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
    E : list
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
    M=ResolMCSVD(A,B)
    [a,b,g]=M
    r=np.sqrt(g+a**2+b**2)
    E=[a,b,r]
    return E

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


def TraceCercle(a,b,r,plot=False):
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
    if plot==True:
        plt.plot(x,y)
        plt.show()
    else :
        return x,y

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
    [[a],[b],[r]]=RegressionLineaireCercle(X,Y)
    x,y=TraceCercle(a,b,r,plot=False)
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
