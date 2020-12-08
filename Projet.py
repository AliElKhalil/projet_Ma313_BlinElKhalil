from Fonctions import *

def ResolMCEN(A,b):
    B = np.dot(A.T,b)
    L = Cholesky(np.dot(A.T,A))
    Lt = L.T #Transposée de L
    Laug = np.column_stack((L,B))
    Y = ResolutionSystTriInf(Laug)
    Ltaug = np.column_stack((Lt,Y))
    X = ResolutionSystTriSup(Ltaug)
    return X

def DecompositionGSr(A):
    n,m=A.shape
    if n !=m :
        raise Exception('Matrice non carrée')

    Q=np.zeros((n,n))
    R=np.zeros((n,n))
    for j in range(n):
        for i in range(j):
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
