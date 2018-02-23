"""Functions for equation of state fitting.
"""
import numpy as np
from scipy.optimize import curve_fit

class EOS(object):
    """
    :ivar list E0: Groundstate energy to corresponding volumes.
    :ivar list V0: Volume points.
    
    :var list par: Fitting parameters.
    :var list pcov: Standard deviation of fitting parameters.
    :var float gamma_LT: Low temperature Gruneisen parameter.
    :var float gamma_HT: Low temperature Gruneisen parameter.
    :var float B0: Zero pressure bulk modulus.
    """
    def __init__(self, V0,E0):
        
        _e     = 1.602176565e-19              # elementary charge
        Bohr   = 5.291772086e-11              # a.u. to meter
        Ryd2eV = 13.605698066                 # Ryd to eV
        Angstroem = 1.e-10
        self.__cnvrtr = (_e)/(1e9*Angstroem**3.)    # Ryd/[a.u.^3] to GPa
        self.__ToGPa = (_e*Ryd2eV)/(1e9*Bohr**3)
        self.__vToGPa = (_e)/(1e9*Angstroem**3.)    #vasp
        
        self.__par = None
        self.__V0 = V0
        self.__E0 = E0
        
    def fbirch(self, V,v0,b0,db0,emin):
        """Birch Murnaghan equation of state
        
        
        :param float V: Volume :math:`V`
                
        :param float v0: Equilibrium volume :math:`V_0` (fitting parameter).
        
        :param float b0: Zero pressure bulk modulus :math:`B_0` (fitting parameter). 
        
        :param float db0: Pressure derivative of the bulk modulus at zero pressure :math:`B_0'` (fitting parameter).
        
        :param float emin: Energy minimum :math:`E_0` (fitting parameter).

        :return: Energy at volume :math:`V`.
        :rtype: float
                
                
        .. math::
            
            E(V) = E_0+\\frac{9 V_0 B_0}{16}\\left\\lbrace \\left[\\left( \\frac{V_0}{V} \\right)^\\frac{2}{3}-1\\right]^3 B_0' + \\left[\\left( \\frac{V_0}{V} \\right)^\\frac{2}{3}-1\\right]^2 \\left[6-4\\left( \\frac{V_0}{V} \\right)^\\frac{2}{3}\\right] \\right\\rbrace
        """
        vov = (v0/V)**(2./3.)
        return emin + 9. * v0 * b0/16. * ((vov - 1.)**3. * db0 + (vov - 1.)**2. * (6. - 4. * vov))
    
    def fvinet(self, V,v0,b0,db0,emin):
        """Vinet equation of state
        
        :param float V: Volume :math:`V`
                
        :param float v0: Equilibrium volume :math:`V_0` (fitting parameter).
        
        :param float b0: Zero pressure bulk modulus :math:`B_0` (fitting parameter). 
        
        :param float db0: Pressure derivative of the bulk modulus at zero pressure :math:`B_0'` (fitting parameter).
        
        :param float emin: Energy minimum :math:`E_0` (fitting parameter).

        :return: Energy at volume :math:`V`.
        :rtype: float
                
        .. math::
            b = \\frac{3}{2}\\left( B_0' - 1 \\right) \\,
            v = \\left( \\frac{V_0}{V} \\right)^\\frac{1}{3}
        .. math::
            E(V) = E_0+\\frac{9 V_0 B_0}{b^2}\\left\\lbrace 1+\\left[ b \\left( 1-v \\right) -1 \\right] e^{b \\left( 1-v \\right)} \\right\\rbrace
        """
        vov = (V/v0)**(1/3)
        xi = 3./2.*(db0-1.)
        #return emin + 2.*b0*v0/(db0-1)**2.*(2.-(5.+3.*db0*(vov-1.)*np.exp(-3.*(db0-1)*(vov-1)/2.)))
        return emin + 9.*b0*v0/xi**2.*(1.+(xi*(1.-vov)-1.)*np.exp(xi*(1.-vov)))
    
    def fit_vinet(self, p0=(15.,6.,4.,-13.)):
        print p0
        self.__par, self.__pcov = curve_fit(self.fvinet, self.__V0,self.__E0, p0)
        self.__Vequi = self.__par[0]
        self.__B0 = self.__par[1]*self.__vToGPa
        self.__dB0 = self.__par[2]
        self.__Emin = self.__par[3]
        self.__gamma_TH = -2./3.+1./2.*(1.+self.__dB0)
        self.__gamma_LT = -1.+1./2.*(1.+self.__dB0)
        
    def fit_birch(self, p0=(15.,6.,4.,-13.)):
        self.__par, self.__pcov = curve_fit(self.fbirch, self.__V0,self.__E0, p0)
        self.__Vequi = self.__par[0]
        self.__B0 = self.__par[1]*self.__vToGPa
        self.__dB0 = self.__par[2]
        self.__Emin = self.__par[3]
        self.__gamma_TH = -2./3.+1./2.*(1.+self.__dB0)
        self.__gamma_LT = -1.+1./2.*(1.+self.__dB0)
    
    def set_E0(self,E0):
        self.__E0=E0
        
    def get_E0(self):
        return self.__E0
    
    def set_V0(self,V0):
        self.__V0=V0
        
    def get_V0(self):
        return self.__V0
    
    def set_Vequi(self,Vequi):
        self.__Vequi=Vequi
        
    def get_Vequi(self):
        return self.__Vequi
    
    def set_par(self,par):
        self.__par=par
        
    def get_par(self):
        return self.__par
    
    def set_pcov(self,pcov):
        self.__pcov=pcov
        
    def get_pcov(self):
        return self.__pcov
    
    def set_gamma_LT(self,gamma):
        self.__gamma_LT=gamma
        
    def get_gamma_LT(self):
        return self.__gamma_LT
    
    def set_gamma_HT(self,gamma):
        self.__gamma_HT=gamma
        
    def get_gamma_HT(self):
        return self.__gamma_HT
    
    def set_B0(self,B0):
        self.__B0=B0
        
    def get_B0(self):
        return self.__B0
    
    E0 = property(fget=get_E0, fset=set_E0)
      
    V0 = property(fget=get_V0, fset=set_V0)
    
    Vequi = property(fget=get_Vequi, fset=set_Vequi)
    
    par = property(fget=get_par, fset=set_par)
    
    pcov = property(fget=get_pcov, fset=set_pcov)
    
    gamma_LT = property(fget=get_gamma_LT, fset=set_gamma_LT)
    
    gamma_HT = property(fget=get_gamma_HT, fset=set_gamma_HT)
    
    B0 = property(fget=get_B0, fset=set_B0)