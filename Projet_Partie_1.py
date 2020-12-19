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
        for i in range(0,j):
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

#def ResolMCQR(A,b):
    #Qr,Rr = DecompositionGSr(A)
    #Taug = np.column_stack((Rr,np.dot(Qr.T,b)))
    #X = ResolutionSystTriSup(Taug)  # RrX = Qrtb
    #Resreturn X

def ResolMCNP(A,b):
    return np.linalg.lstsq(A,b)[0]


def DecompositionGSGenerale(A):
    """
    

    Parameters
    ----------
    A : Matrice 
        Matrice que l'on souhaite décomposer sous la forme QR. 
        De format (n,p), on peut avoir n=p, n>p et n<p.

    Raises
    ------
    Exception
        Décomposition QR impossible si on n'a pas Ker(A)={0}.

    Returns
    -------
    Q : Matrice 
        Q, de taille (n,p), vérifie Qt.Q=Ip, matrice identitée de taille p.
        Attention : Q.Qt n'est pas nécessairement égale à In.
    R : Matrice
        Matrice triangulaire supérieure de taille p.

    """
    
    n,p=np.shape(A)
    Q=np.zeros((n,p))
    R=np.zeros((p,p))
    for j in range (0,p):
        for i in range (0,j):
            R[i,j]=np.vdot(A[:,j],Q[:,i])
        w=A[:,j]
        for k in range(j):
            w=w-R[k,j]*Q[:,k]
        norme=np.linalg.norm(w) 
        if norme ==0:
            raise Exception('décomposition QR impossible : Ker(A)!={0}')
        R[j,j]=norme
        Q[:,j]=(1/norme)*w
    return Q,R

def ResolMCQR(A,b):
    Q,R=DecompositionGSGenerale(A)
    Qt=Q.T
    x=ResolTriSup(R,np.dot(Qt,b))
    return x
