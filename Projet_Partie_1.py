from Fonctions import *
import numpy as np


def ResolMCEN(A,b):
    B = np.dot(A.T,b)
    L = Cholesky(np.dot(A.T,A))
    Lt = L.T #Transposée de L
    #Laug = np.column_stack((L,B))
    #Y = ResolutionSystTriInf(Laug)
    #Ltaug = np.column_stack((Lt,Y))        
    #Ltaug = np.column_stack((Lt,Y))
    #X = ResolutionSystTriSup(Ltaug)
    Y=ResolTriInf(L,B)
    X=ResolTriSup(Lt,Y)
    return X

def DecompositionGSr(A):
    n,p=A.shape
    if n !=p :
        print ('Décomposition GS réduite')
    if n < p :
        raise Exception ("n doit être >= p")

    Q = np.zeros((n,p))
    R = np.zeros((n,p))
    for j in range(p):
        for i in range(j,p):
            R[i,j]=Q[:,i]@A[:,j]
        w=A[:,j]
        for k in range(j):
            w=w-R[k,j]*Q[:,k]
        norme=np.linalg.norm(w)
        if norme ==0:
            raise Exception('Matrice non inversible')
        R[j,j]=norme
        Q[:,j]=w/norme
    return Q,R

def ResolMCQR(A,b):
    Qr,Rr = DecompositionGSr(A)
    Taug = np.column_stack((Rr,np.dot(Qr.T,b)))
    X = ResolutionSystTriSup(Taug)  # RrX = Qrtb
    return X

def ResolMCNP(A,b):
    X = "OUIIII"
    return X


def DecompositionGSGenerale(A):
    n,p=np.shape(A)
    Q=np.zeros((n,p))
    R=np.zeros((p,p))
    R[0,0]=np.linalg.norm(A[:,0])    
    Q[:,0]=(1/R[0,0])*A[:,0]
    for j in range (1,p):
        for i in range (1,j):
            R[i,j]=np.vdot(A[:,j],Q[:,i])
        w=A[:,j]
        for k in range(j):
            w=w-R[k,j]*Q[:,k]
        norme=np.linalg.norm(w)
        if norme ==0:
            raise Exception('décomposition QR impossible : élément diagonle de R nul')
        R[j,j]=norme
        Q[:,j]=w/norme
    return Q,R


