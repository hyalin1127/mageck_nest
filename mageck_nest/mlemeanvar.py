'''
Modeling mean and variance
'''
import os
import re
import scipy
from scipy.stats import nbinom,linregress
import random
import math
import numpy as np
import numpy.linalg as linalg
import logging
from scipy.optimize import curve_fit
import statsmodels.api as sm
import statsmodels.formula.api as smf
import operator
import matplotlib.pyplot as plt

# debug
try:
  from IPython.core.debugger import Tracer
except:
  pass

from mageck_nest.mleclassdef import *

def func(x,a,b):
    return (a/x)+b

class MeanVarModel:
    '''
    The Mean and Dispersion model class
    '''
    # Ror linear regression
    #lm_intercept=0.0
    #lm_coeff=0.0

    # For generalized linear model
    glm_a0=0.0
    glm_a1=0.0
    '''
    def model_mean_disp_by_lm(self,allgenedict):
        '
        #Modeling the mean and dispersion by linear regression

        list_k=[]
        list_dispersion=[]
        for (gid,gsk) in allgenedict.iteritems():
            nsg=len(gsk.nb_count[0])
            nsample=len(gsk.nb_count)
            if len(gsk.sgrna_kvalue)>0:
                if gsk.MAP_sgrna_dispersion_estimate!=None:
                    sg_k=[x[0] for x in gsk.sgrna_kvalue.tolist()]
                    sg_dispersion=gsk.MAP_sgrna_dispersion_estimate
                    if len(sg_k)>=nsg*nsample:
                        list_k+=sg_k[:(nsg*nsample)]
                        list_dispersion+=sg_dispersion[:(nsg*nsample)]


        k_log=np.log(list_k)
        dispersion_log=np.log(list_dispersion)
        # remove those with too low variance
        k_log2=np.array([k_log[i] for i in range(len(dispersion_log)) if dispersion_log[i]>(-1)])
        dispersion_log2=np.array([dispersion_log[i] for i in range(len(dispersion_log)) if dispersion_log[i]>(-1)])
        if len(k_log2)>20:
            (slope,intercept,r_value,p_value,std_err)=linregress(k_log2,dispersion_log2)
        else:
            (slope,intercept,r_value,p_value,std_err)=linregress(k_log,dispersion_log)
        self.lm_intercept=intercept
        self.lm_coeff=slope

        logging.info('Linear regression: y='+str(slope)+'x+'+str(intercept))
        if np.isnan(slope) or np.isnan(intercept):
            logging.error('Nan values for linear regression')

    def get_lm_dispersion(self,klist,returnalpha=False):

        #Return the fitted values of variance.
        #If returnalpha=True, return the alpha value (var=mean+alpha*mean^2)

        kls=(klist)
        k_log=np.log(kls)
        dispvalue=k_log*self.lm_coeff+self.lm_intercept
        dispvalue=np.exp(dispvalue)
        return(dispvalue)
    '''
    def model_mean_disp_by_glm(self,allgenedict,output_prefix,size_f):
        '''
        Fitting the mean variance to the theoreotical curve
        '''
        list_k=[]
        list_dispersion=[]

        inverse_size_f=[1/i for i in size_f]
        for (gid,gsk) in allgenedict.items():
            nallsample=gsk.nb_count.shape[0]
            n=gsk.nb_count.shape[1]
            if gsk.MAP_sgrna_dispersion_estimate!=None:
                sg_dispersion=gsk.MAP_sgrna_dispersion_estimate
                sg_k=[x[0] for x in gsk.sgrna_kvalue.tolist()]
                list_dispersion+=sg_dispersion[:len(gsk.sgrnaid)]
                for i in range(len(gsk.sgrnaid)):
                    normalized_k_mean=np.mean(np.multiply(inverse_size_f,[sg_k[i+k*n] for k in range(nallsample)]))
                    #normalized_k_mean=np.mean(np.multiply(inverse_size_f,sg_k[i*nallsample:(i+1)*nallsample]))
                    list_k.append(normalized_k_mean)

        combined_list=[[list_k[i],list_dispersion[i]] for i in range(len(list_k)) if math.isnan(list_k[i])==False and list_dispersion[i]<10 and math.isnan(list_dispersion[i])==False]
        xdata=np.array([i[0] for i in combined_list])
        ydata=np.array([i[1] for i in combined_list])
        popt, pcov = curve_fit(func, xdata, ydata)

        self.glm_a0=popt[0]
        self.glm_a1=popt[1]

        xdemo = np.linspace(0, 2000, 500)
        ydemo = func(xdemo,popt[0],popt[1])

    def get_glm_dispersion(self,klist,returnalpha=False):
        kls=np.array(klist)
        dispvalue=func(kls,self.glm_a0,self.glm_a1)
        return(dispvalue)
