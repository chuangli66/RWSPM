"""
simulate data
  for 1 dataset
"""

import os
import numpy as np
from scipy.stats import norm, binom
import math
from scipy.io import savemat
import time
import pyreadr
import pandas as pd
import os.path


k = 7

n = 200
a0 = 0.25
sigma = 3


M1 = 150
M2 = 100
M = M1*M2

""" affected region """
c = 50 # x0 = y0 = c: true center of affected region
r0 = 5 # radius of affected region
sigma_r = 1 # standard error of r0

for m in range(1,101):
    data_folder = ".data/result/Setting" + str(k) + "/sim_" + str(m)
    if not os.path.exists(data_folder):
       os.mkdir(data_folder)
    
    start_2 = time.time()
    print(m)
    
    """ epsilon """
    # np.random.seed(seed=m)
    xi_10 = norm.rvs(size=n)
    xi_1 = np.sign(xi_10) * (np.abs(xi_10) - np.floor(np.abs(xi_10)))
    xi_2 = np.sqrt(1-xi_1**2)*(binom.rvs(n=1, p=0.5, size=n)*2 - 1)
    
    
    epsilon = np.zeros(shape = (n,M))
    
    
    for i in np.arange(n):
    	base = np.zeros(shape = (M1,M2))
    	for x in np.arange(M1):
    		for y in np.arange(M2):
    			psi_1 = 2**0.5 * np.sin(2*math.pi*(x+1)/M1)
    			psi_2 = 2**0.5 * np.cos(2*math.pi*(y+1)/M2)
    			base[x,y] = xi_1[i]*psi_1 + xi_2[i]*psi_2
    	epsilon[i,:] = base.reshape((1,M))
    
    # np.min(epsilon): -1.9997
    # np.max(epsilon):  1.9997
    # np.var(epsilon): 1
    """ Z """
    # np.random.seed(seed=m)
    Z1 = np.random.uniform(50, 100, size = n)
    Z2 = binom.rvs(n=1, p=0.5, size=n)
    
    """ X """
    " X ~ Bin(2, 0.2), for n=200 "
    # np.random.seed(seed=m)
    X = binom.rvs(n=2, p=0.2, size=(n,1))
    
    """ Y """
    prod = np.zeros(shape = (n,M)) # 200 x 15000
    # np.random.seed(seed=m)

    c1 = c + (sigma*norm.rvs(size=n)).astype(int)
    c2 = c + (sigma*norm.rvs(size=n)).astype(int)
    r =  np.multiply(r0 + sigma_r*norm.rvs(size=n), X.flatten()) + Z1/100 + Z2
    
    time0 = time.time()

    gp1 = np.arange((c - 10),(c + 10))
    gp2 = np.arange((c - 10),(c + 10))
    for i in np.arange(n):
      idx = np.zeros(shape = (M1,M2))
      for x in gp1:
        for y in gp2:
          if (x-c1[i])**2 + (y-c2[i])**2 <= r[i]**2:
            # beta = 0.5*np.exp(-0.1*np.sqrt(((x-c1[i])**2) + ((y-c2[i])**2)))
            beta = 0.8*np.exp(-0.1*np.sqrt(((x-c1[i])**2) + ((y-c2[i])**2)))
            idx[x,y] = (X[i][0] * beta)**2 + 0.01*Z1[i] + 0.25*Z2[i]
      prod[i,:] = idx.reshape((1,M))
    time1 = time.time() - time0
    
    Y = prod + epsilon
    # min(y0) = min(y1) = min(y2) = -1.99973
    # max(y0) = 1.999732; max(y1) = 2.27; max(y2) =  2.77
    
    """ save """

    file_name = '%s/x.RData' % data_folder
    df = pd.DataFrame(X)
    pyreadr.write_rdata(file_name, df, df_name="x")
    
    file_name = '%s/epsilon.RData' % data_folder
    df = pd.DataFrame(epsilon)
    pyreadr.write_rdata(file_name, df, df_name="epsilon")
    
    file_name = '%s/psrod.RData' % data_folder
    df = pd.DataFrame(prod)
    pyreadr.write_rdata(file_name, df, df_name="prod")
    
    file_name = '%s/y.RData' % data_folder
    df = pd.DataFrame(Y)
    pyreadr.write_rdata(file_name, df, df_name="y")
    
    file_name = '%s/xi1.RData' % data_folder
    df = pd.DataFrame(xi_1)
    pyreadr.write_rdata(file_name, df, df_name="xi_1")
    
    file_name = '%s/xi2.RData' % data_folder
    df = pd.DataFrame(xi_2)
    pyreadr.write_rdata(file_name, df, df_name="xi_2")
    
    file_name = '%s/Z1.RData' % data_folder
    df = pd.DataFrame(Z1)
    pyreadr.write_rdata(file_name, df, df_name="Z1")
    
    file_name = '%s/Z2.RData' % data_folder
    df = pd.DataFrame(Z2)
    pyreadr.write_rdata(file_name, df, df_name="Z2")
    
    end_2 = time.time()
    print("Elapsed time is ", round((end_2 - start_2)), " seconds.") # 29 seconds




