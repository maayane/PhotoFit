
"""*******************************************************
Classes and functions for model fitting
******************************************************"""
__author__ = 'maayanesoumagnac'


import numpy as np
from scipy import optimize


class objective_no_cov(object):
    def __init__(self,interpolated_model,data): #the inversion is performed outside of the class to save time in loop
        self.interpolated_model=interpolated_model#array of X,Y model
        self.data=data#array of X,Y data
    def Res(self):
#        print 'I am calculating the Res'
        #print np.ndim(self.data)
        #pdb.set_trace()
        if np.ndim(self.data)>1:
            res=self.interpolated_model[:]-self.data[:,1]#res- a numpy array[model(interpolated)-data]
        else:
            res = self.interpolated_model[:] - self.data[:]  # res- a numpy array[model(interpolated)-data]
        return res
    def chi_square_value(self):
        chi2=np.dot(self.Res().transpose(),self.Res())
        return chi2

class objective_with_cov(object):
    def __init__(self,interpolated_model,data,inverted_covariance): #the inversion is performed outside of the class to save time in loop
        self.interpolated_model=interpolated_model#Y model
        self.data=data#array of X,Y data
        self.inverted_covariance=inverted_covariance #array of covariance
    def Res(self):
        if np.ndim(self.data)>1:
            #print('np.ndim(self.data)>1')
            res=self.interpolated_model[:]-self.data[:,1]#res- a numpy array[model(interpolated)-data]
            #print('res from np.ndim(self.data)>1 is',res)
            #print('data from np.ndim(self.data)>1 is',self.data[:,1])
        else:
            res = self.interpolated_model[:] - self.data[:]
            #print('res from else is',res)
            #print('data from else is',self.data[:])

        return res
    def chi_square_value(self):
        #print('np.shape(self.Res().transpose()) is',np.shape(self.Res().transpose()))
        #print('np.shape(self.inverted_covariance) is',np.shape(self.inverted_covariance))
        #print('np.shape(np.dot(self.inverted_covariance,self.Res())) is',np.shape(np.dot(self.inverted_covariance,self.Res())))
        chi2=np.dot(self.Res().transpose(),np.dot(self.inverted_covariance,self.Res()))
        #print('chi2 is',chi2)
        return chi2

class Parameter:
    def __init__(self, value):
            self.value = value

    def set(self, value):
            self.value = value

    def __call__(self):
            return self.value

def fit(function, parameters, y, x = None):
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return y - function(x)
    if x is None: x = np.arange(y.shape[0])
    p = [param() for param in parameters]
    #print('x is',x)
    #print('function(x) is',function(x))
    #print('y-function(x) is',y-function(x))
    #print('f(p) is',f(p))
    #print('p is',p)
    #print('optimize.leastsq(f,p) is',optimize.leastsq(f,p))
    return optimize.leastsq(f, p)
    #the way to run it is:
        # giving initial parameters
        #mu = Parameter(7)
        #sigma = Parameter(3)
        #height = Parameter(5)
        # define your function:
        # def f(x): return height() * np.exp(-((x-mu())/sigma())**2)
        # print fit(f, [mu, sigma, height], data)




