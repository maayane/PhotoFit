#! //anaconda/bin/python

"""*******************************************************
This code contains classes usefull for fitting
******************************************************"""
#print __doc__

__author__ = 'maayanesoumagnac'



import numpy as np
import scipy as sp
from scipy import interpolate



class objective_with_uncertainties_notint(object):#the case where the model and the data do not have the same x-axix
        def __init__(self, model, data, sigmas):
            self.model = model # a 2-n array
            self.data = data # a 2-m array
            self.sigmas = sigmas # a m-long vector, with the uncertainties sigmas
            #self.intmodel=intmodel # fasten the code. If not None, an 2-m array with model interpolated on the data x axis
        def Res(self):
            model=self.model
            interpolate_model=np.zeros((np.shape(self.data)))
            Intmodel=sp.interpolate.interp1d(self.model[:, 0], self.model[:, 1])
            interpolate_model[:,0]=self.data[:,0]
            interpolate_model[:,1]=Intmodel(self.data[:,0])
            Res = np.abs(interpolate_model[:,1] - self.data[:,1])# in the absence of cov, the Res will be elevated to square anyway
            return Res
        def chi_square_value(self):
            invcov=np.diag(np.power(self.sigmas,-2))
            chi2 = np.dot(self.Res().transpose(), np.dot(invcov, self.Res()))
            return chi2


class objective_with_uncertainties(object):#the case where the model and the data do not have the same x-axix
    def __init__(self, model, data, sigmas):
        self.model = model  # a 2-n array
        self.data = data  # a 2-m array
        self.sigmas = sigmas  # a m-long vector, with the uncertainties sigmas
        #self.intmodel = intmodel  # fasten the code. If not None, an 2-m array with model interpolated on the data x axis
    def Res(self):
        interpolate_model = self.model
        #print interpolate_model[:,1]
        #print self.data[:,1]
        #pdb.set_trace()
        Res = np.abs(interpolate_model[:, 1] - self.data[:, 1])# in the absence of cov, the Res will be elevated to square anyway
        #print 'Res is',Res
        return Res
    def chi_square_value(self):
        invcov = np.diag(np.power(self.sigmas, -2))
        chi2 = np.dot(self.Res().transpose(), np.dot(invcov, self.Res()))
        #print 'chi2 is',chi2
        return chi2