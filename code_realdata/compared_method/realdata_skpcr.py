
import numpy as np

from scipy.spatial.distance import cdist
import scipy.sparse.linalg
from scipy.stats import rankdata,norm, binom
import scipy
import time
import os
import pyreadr

def sKPCR_prepare_Laplacian(mat):
    #mask:   the nxp 2D matrix (n=200, p=21)
    #output: the laplacian matrix
    p = mat.shape[1]
    dims=np.vstack((range(p),range(p)))
    dist_map=cdist(dims.T,dims.T,'euclidean')
    tmp=np.where((dist_map<2) & (dist_map>0))
    adj=np.zeros((p,p))
    adj[tmp]=1
    degree = np.diag([sum(x) for x in zip(*adj)])
    laplacian_mat=degree-adj
    return laplacian_mat


def sKPCR_linear(X,laplacian_mat,K):
    #X is n*p, laplacian_mat is p*p sparse, k is a number (5)
    n=X.shape[0]
    Kernel0=np.dot(np.dot(X,laplacian_mat),X.T)
    oneN=np.divide(np.ones((n,n)),np.float64(n))
    Kernel=np.divide((Kernel0-np.dot(oneN,Kernel0)-np.dot(Kernel0,oneN)+np.dot(np.dot(oneN,Kernel0),oneN)),np.float64(n)) # center, /n
    D,U=np.linalg.eig(Kernel)
    D=np.real(D)
    U=np.real(U)
    indx1=np.argsort(-D)
    D=D[indx1]
    U=U[:,indx1]
    U=U[:,0:K]
    return U



def BWAS_regression(X,Y):
    # X=yperms[:,j]
    # Y=pcs
    n,P = Y.shape
    df=np.float32(X.shape[0]-1) #degree of freedom
    ss=np.dot(X.T,X)
    beta=np.dot(np.dot(ss,np.transpose(X)),Y)
    Res=Y-np.dot(X.reshape((n,1)),beta.reshape((1,P)))
    sigma=np.reshape(np.sqrt(np.divide(np.sum(np.square(Res),axis=0),df)),(1,P))
    Tstat=np.divide(np.transpose(beta),np.dot(sigma,np.sqrt(ss)))
    return Tstat


def sKPCR_adptive_regression(pcs,pheno,numperms):
    N,P=pcs.shape
    # yperms: 
    #   the first one is pheno
    #   the 2nd to end are permuted pheno
    
    yperms = np.zeros((N,numperms),dtype='float32')
    yperms[:,0]=pheno
    for j in range(1,numperms):
        yperms[:,j] = pheno[np.random.permutation(N)]
        
    Tstat2=np.zeros((numperms,P),dtype='float32')
    for j in range(0,numperms):
        Tstat2[j,:]=np.square(BWAS_regression(yperms[:,j],pcs))
    
    pPerm0=np.ones((P,1),dtype='float32') 
    T2stats=np.ones((P,),dtype='float32') 
    for j in range(0,P):
        #test statistic (1+numperm)*1 vector
        T0s = np.mean(Tstat2[:,0:(j+1)],1)
        T2stats[j]=T0s[0]
        #pvalue
        pPerm0[j]= np.mean(T0s[0]<=T0s[1:numperms])
        #get ranks of 
        ranks=rankdata(T0s[1:numperms])
        #get empirical p-value using rank
        P0s = (numperms-ranks)/(np.float32(numperms-1))
        if j==0:
            minp0=P0s*1.0
        else:
            minp0[minp0>P0s]=P0s[minp0>P0s]*1.0
    
    pvs=max(1/np.float32(numperms),(np.sum(minp0<=np.min(pPerm0)))/np.float32(numperms-1))
    T2stats_best=T2stats[np.where(np.min(pPerm0)==pPerm0)[0][0]]
    
    return  pvs,T2stats_best

def sKPCR_run_full_analysis(mat, x, K_components,numperms):
    # start_time = time.time()
    laplacian_mat=sKPCR_prepare_Laplacian(mat)
    U=sKPCR_linear(mat,laplacian_mat,K_components)
    # numperms = 99
    pv,T2stat=sKPCR_adptive_regression(U,x,numperms)
    # end_time = time.time()
    # t1 = end_time - start_time
    return pv



### global test ###
x1 = pyreadr.read_r('/.../ADNI_1_3.Rdata')['ADNI_COMBINE'].values[:, 1:]
### left hipp 
# y = pyreadr.read_r('/.../left_res_ori.Rdata')['left_res_ori']
### right hipp 
y = pyreadr.read_r('/data/right_res_ori.Rdata')['right_res_ori']
o = np.zeros((494465, 2))
for l in range(1,494466):
    pvals = sKPCR_run_full_analysis(y, x1[:,l], 5,99)
    if  pvals < 0.02:
        pvals = sKPCR_run_full_analysis(y, x1[:,l], 5,9999)
        if  pvals < 0.002:
            pvals = sKPCR_run_full_analysis(y, x1[:,l], 5,99999)
    end = time.time()
    o[l - 1, 0] = l
    o[l - 1, 1] = pvals
# np.savetxt('/.../skpcr_global_left.txt', o)
np.savetxt('/.../skpcr_global_right.txt', o)

