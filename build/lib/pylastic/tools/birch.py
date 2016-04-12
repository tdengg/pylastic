import numpy as np
from scipy.optimize import curve_fit

class EOS(object):
    def __init__(self, V0,E0):
        self.__par = None
        self.__V0 = V0
        self.__E0 = E0
        
    def fbirch(self, V,v0,b0,db0,emin):
        vov = (v0/V)**(2./3.)
        return emin + 9. * v0 * b0/16. * ((vov - 1.)**3. * db0 + (vov - 1.)**2. * (6. - 4. * vov))

    def fvinet(self, V,v0,b0,db0,emin):
        vov = (V/v0)**(1./3.)
        xi = 3./2.*(db0-1.)
        #return emin + 2.*b0*v0/(db0-1)**2.*(2.-(5.+3.*db0*(vov-1.)*np.exp(-3.*(db0-1)*(vov-1)/2.)))
        return emin + 9.*b0*v0/xi**2.*(1.+(xi*(1.-vov)-1.)*np.exp(xi*(1.-vov)))
    def fit_vinet(self):
        self.__par, self.__pcov = curve_fit(self.fvinet, self.__V0,self.__E0, p0=(15.,6.,4.,-13.))
    def fit_birch(self):
        self.__par, self.__pcov = curve_fit(self.fbirch, self.__V0,self.__E0, p0=(15.,6.,4.,-13.))
    
    def set_E0(self,E0):
        self.__E0=E0
        
    def get_E0(self):
        return self.__E0
    
    def set_V0(self,V0):
        self.__V0=V0
        
    def get_V0(self):
        return self.__V0
    
    def set_par(self,par):
        self.__par=par
        
    def get_par(self):
        return self.__par
    
    E0 = property(fget=get_E0, fset=set_E0)
      
    V0 = property(fget=get_V0, fset=set_V0)
    
    par = property(fget=get_par, fset=set_par)