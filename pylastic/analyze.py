""" Determination of elastic constant and cross validation score.
"""

import numpy as np
import math
from copy import copy
import matplotlib.pyplot as plt

class Energy(object):
    """ Implementation of the energy aproach for determination of elastic constants.
    """
    def __init__(self, strain, energy, V0):
        _e        =  1.602176565e-19              # elementary charge
        Bohr      =  5.291772086e-11              # a.u. to meter
        Ryd2eV    = 13.605698066                  # Ryd to eV
        Angstroem =  1.e-10                       # Angstroem to meter
        self.__vToGPa = (_e)/(1e9*Angstroem**3.)
        
        self.__V0 = V0
        self.__strain = strain
        self.__energy = energy
        self.__Cij2nd  = {}
        
        
    def set_2nd(self, fitorder):
        """Fit energy strain curve and evaluate 2nd derivative in order to get 2nd order elastic constants."""
        self.__CONV = self.__vToGPa #* math.factorial(2)*2.
        strain = copy(self.__strain)
        energy = copy(self.__energy)
        while (len(strain) > fitorder): 
            
            emax  = max(strain)
            emin  = min(strain)
            emax  = max(abs(emin),abs(emax))
            coeffs= np.polyfit(strain, energy, fitorder)
            
            self.__Cij2nd[str(emax)]  = coeffs[fitorder-2]*self.__CONV/self.__V0         # in GPa units 
            
            if (abs(strain[0]+emax) < 1.e-7):
                strain.pop(0); energy.pop(0)
            if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                strain.pop()
                energy.pop()
        
    def get_2nd(self):
        return self.__Cij2nd
    
    
    def set_3rd(self, fitorder):
        """Evaluate 3rd order elastic constants from energy strain curves"""
        self.__CONV = self.__vToGPa * math.factorial(3)*2.
        strain = copy(self.__strain)
        energy = copy(self.__energy)
        while (len(strain) > fitorder): 
            
            emax  = max(strain)
            emin  = min(strain)
            emax  = max(abs(emin),abs(emax))
            coeffs= np.polyfit(strain, energy, fitorder)
            self.__Cij3rd  = coeffs[fitorder-3]*self.__CONV/self.__V0 * 0.001 # in TPa units
            
            if (abs(strain[0]+emax) < 1.e-7):
                strain.pop(0); energy.pop(0)
            if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                strain.pop()
                energy.pop()
        
    def get_3rd(self):
        return self.__Cij3rd
    
    
    def set_cvs(self, fitorder):
        """Evaluation of cross validation score for fitorder n."""
        strain = copy(self.__strain)
        energy = copy(self.__energy)
        
        while (len(strain) > fitorder+1):
            emax = max(strain)
            emin = min(strain)
            emax = max(abs(emin),abs(emax))
    
            S = 0
            for k in range(len(strain)):
                Y      = energy[k]
                etatmp = []
                enetmp = []
    
                for l in range(len(strain)):
                    if (l==k): pass
                    else:            
                        etatmp.append(strain[l])
                        enetmp.append(energy[l])
    
                Yfit = np.polyval(np.polyfit(etatmp,enetmp, fitorder), strain[k])
                S    = S + (Yfit-Y)**2
            
            self.__CV = np.sqrt(S/len(strain))
            
    
            if (abs(strain[0]+emax) < 1.e-7):
                strain.pop(0)
                energy.pop(0)
            if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                strain.pop()
                energy.pop()
                
                
    def get_cvs(self):
        return self.__CV
    
    def plot_energy(self):
        """Return matplotlib axis instance for energy-strain curve.  """
        
        ax = plt.plot(self.__strain, self.__energy)
        return ax
        

        
#class Strain(object):
#    def __init__(self):
        