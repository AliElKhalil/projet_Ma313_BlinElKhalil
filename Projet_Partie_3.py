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
    X=[X]
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
        Coordonnees en ordonee des points a etudier.

    Returns
    -------
    None.

    """
    [[a],[b],[g]]=RegressionLineaireCercle(X,Y)
    R=np.sqrt(g+a**2+b**2)
    theta=np.linspace(0,2*np.pi,1000)
    x=R*np.cos(theta)+a
    y=R*np.sin(theta)+b
    plt.plot(x,y)
    plt.scatter(X,Y)
    plt.show()