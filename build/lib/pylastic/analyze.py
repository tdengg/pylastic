""" Determination of elastic constant and cross validation score.
"""

import numpy as np
import math
from copy import copy
try:
    import matplotlib.pyplot as plt
    mpl=True
except:
    mpl=False
    
class Energy(object):
    """ Implementation of the energy aproach for determination of elastic constants.
    
    Parameters
    ----------
    strain : list
        List of energy values (float).
    strain : list
        List of strain values (float).
    V0 : float 
        Equilibrium volume of parent structure.
    """
    def __init__(self, strain=None, energy=None, V0=None, code='vasp'):
        _e        =  1.602176565e-19              # elementary charge
        Bohr      =  5.291772086e-11              # a.u. to meter
        Ryd2eV    = 13.605698066                  # Ryd to eV
        Angstroem =  1.e-10                       # Angstroem to meter
        self.__vToGPa = (_e)/(1e9*Angstroem**3.)
        self.__ToGPa  = (_e*Ryd2eV)/(1e9*Bohr**3.)
        
        self.__V0 = V0
        self.__cod = code
        self.__strain = strain
        self.__energy = energy
        self.__Cij2nd  = {}
        self.__CV = []
        self.__coeffs = {}
        self.__r = {}
        
    def set_strain(self, strain):
        self.__strain=strain
    
    def get_strain(self):
        return self.__strain
    
    def set_energy(self, energy):
        self.__energy = energy
    
    def get_energy(self):
        return self.__energy
    
    def set_V0(self,V0):
        self.__V0 = V0
        
    def get_V0(self):
        return self.__V0
        
    def set_2nd(self, fitorder):
        """Fit energy vs. strain curve and evaluate 2nd derivatives in order to get 2nd order elastic constants.
        
        Parameters
        ----------
        fitorder : integer 
            Order of polynomial energy-strain fit.
        """
        
        self.search_for_failed()
        #self.__CONV = self.__vToGPa #* math.factorial(2)*2.
        if self.__cod == 'vasp': self.__CONV = self.__vToGPa * math.factorial(2)*2.
        if self.__cod == 'wien': self.__CONV = self.__ToGPa * math.factorial(2)*1.
        if self.__cod == 'espresso': self.__CONV = self.__ToGPa * math.factorial(2)*1.
        if self.__cod == 'exciting': self.__CONV = self.__ToGPa * math.factorial(2)*2.
        if self.__cod == 'emto': self.__CONV = self.__ToGPa * math.factorial(2)*2.
        strain = copy(self.__strain)
        energy = copy(self.__energy)
        
        while (len(strain) > fitorder): 
            
            emax  = max(strain)
            emin  = min(strain)
            emax  = max(abs(emin),abs(emax))
            coeffs= np.polyfit(strain, energy, fitorder)
            
            self.__Cij2nd[str(emax)]  = coeffs[fitorder-2]*self.__CONV/self.__V0         # in GPa units 
            print self.__V0
            """  Calculate RMS:  """
            deltasq = 0
            deltas = []
            
            poly = np.poly1d(coeffs)
            
            for i in range(len(energy)):
                delta = energy[i] - poly(strain[i])
                deltas.append(delta)
                deltasq += (delta)**2.0
                
            self.__r[(emax,fitorder)] = np.sqrt(deltasq/len(strain))
            
            self.__coeffs[(round(emax,4),fitorder)] = (coeffs)
            
            if (abs(strain[0]+emax) < 1.e-7):
                strain.pop(0); energy.pop(0)
            if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                strain.pop()
                energy.pop()
                
    def get_2nd(self):
        return self.__Cij2nd
    
    def get_r(self):
        return self.__r
    
    def set_3rd(self, fitorder):
        """Evaluate 3rd order elastic constants from energy strain curves.
        
        Parameters
        ----------
        fitorder : integer 
            Order of polynomial energy-strain fit.
        """
        self.search_for_failed()
        
        if self.__cod == 'vasp': self.__CONV = self.__vToGPa * math.factorial(3)*2.
        if self.__cod == 'wien': self.__CONV = self.__ToGPa * math.factorial(3)*1.
        if self.__cod == 'espresso': self.__CONV = self.__ToGPa * math.factorial(3)*1.
        if self.__cod == 'exciting': self.__CONV = self.__ToGPa * math.factorial(3)*2.
        if self.__cod == 'emto': self.__CONV = self.__ToGPa * math.factorial(3)*2.
        strain = copy(self.__strain)
        energy = copy(self.__energy)
        while (len(strain) > fitorder): 
            
            emax  = max(strain)
            emin  = min(strain)
            emax  = max(abs(emin),abs(emax))
            self.__coeffs= np.polyfit(strain, energy, fitorder)
            self.__Cij3rd  = self.__coeffs[fitorder-3]*self.__CONV/self.__V0 * 0.001 # in TPa units
            
            """  Calculate RMS:  """
            deltasq = 0
            deltas = []
            
            poly = np.poly1d(self.__coeffs)
            
            for i in range(len(energy)):
                delta = energy[i] - poly(strain[i])
                deltas.append(delta)
                deltasq += (delta)**2.0
                
            self.rms.append(np.sqrt(deltasq/len(strain)))
            
            if (abs(strain[0]+emax) < 1.e-7):
                strain.pop(0); energy.pop(0)
            if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                strain.pop()
                energy.pop()
                
        
    def get_3rd(self):
        return self.__Cij3rd
    
    
    def set_cvs(self, fitorder):
        """Evaluation of cross validation score.
        
        Parameters
        ----------
        fitorder : integer 
            Order of polynomial energy-strain fit.
        """
        self.search_for_failed()
        
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
            
            self.__CV.append((np.sqrt(S/len(strain)),emax,fitorder))
            
    
            if (abs(strain[0]+emax) < 1.e-7):
                strain.pop(0)
                energy.pop(0)
            if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                strain.pop()
                energy.pop()
                
                
    def get_cvs(self):
        return self.__CV

    def search_for_failed(self, mod='pass'):
        """Handle failed calculations"""
        i=0
        while i < len(self.__strain):
            
            if self.__energy[i] == 0:
                print "WARNING: Calculation failed on %s ...... trying to fit anyway."%str(self.__strain[i])
                self.__energy.pop(i)
                self.__strain.pop(i)
            i+=1
            
    energy = property( fget = get_energy        , fset = set_energy    )
    strain = property( fget = get_strain        , fset = set_strain    )
    V0 = property( fget = get_V0        , fset = set_V0    )
                
        
class Stress():
    """ Implementation of the stress aproach for determination of elastic constants.
    
    Parameters
    ----------
    strain : list
        List of energy values (float).
    strain : list
        List of strain values (float).
    V0 : float 
        Equilibrium volume of parent structure.
    """
    def __init__(self, strain=None, stress=None, V0=None, code='vasp'):
        _e        =  1.602176565e-19              # elementary charge
        Bohr      =  5.291772086e-11              # a.u. to meter
        Ryd2eV    = 13.605698066                  # Ryd to eV
        Angstroem =  1.e-10                       # Angstroem to meter
        self.__vToGPa = (_e)/(1e9*Angstroem**3.)
        
        self.__V0 = V0
        self.__cod = code
        self.__strain = strain
        self.__stress = stress
        self.__Cij2nd  = {}
        self.__CV = []
        self.__coeffs = {}
        self.__r = {}
        
    def set_strain(self, strain):
        self.__strain=strain
    
    def get_strain(self):
        return self.__strain
    
    def set_stress(self, stress):
        self.__stress = stress
    
    def get_stress(self):
        return self.__stress
    
    def set_V0(self,V0):
        self.__V0 = V0
        
    def get_V0(self):
        return self.__V0
        
    def set_2nd(self, fitorder):
        """Fit stress strain curve and evaluate 2nd derivative in order to get 2nd order elastic constants.
        
        Parameters
        ----------
        fitorder : integer 
            Order of polynomial stress-strain fit.
        """
        
        
        
        
        self.search_for_failed()
        self.__CONV = self.__vToGPa #* math.factorial(2)*2.
        #strain = copy(self.__strain)
        #stress_all = copy(self.__stress)
        
        LSi_dic = {1:'LS1',2:'LS2',3:'LS3',4:'LS4',5:'LS5',6:'LS6'}
        sigma = {}
        sigma1 = []
        sigma2 = []
        sigma3 = []
        sigma4 = []
        sigma5 = []
        sigma6 = []
        
        for l in range(1,7):

            
            if l==1: s=(0,0)
            elif l==2: s=(1,1)  
            elif l==3: s=(2,2) 
            elif l==4: s=(1,2) 
            elif l==5: s=(0,2) 
            elif l==6: s=(0,1) 
            stress_ii = [st[s[0]][s[1]] for st in self.__stress]
            
            strain = copy(self.__strain)
            stress = copy(stress_ii)
            print strain, stress, s
            etacalc = []
            #--- first derivative coefficient calculation -----------------------------------------
            
            while (len(strain) > fitorder):
                emax=max(strain)
                emin=min(strain)
                emax=max(abs(emin),abs(emax))
                
                coeffs = np.polyfit(strain, stress, fitorder)
                
                #self.__Cij2nd[str(emax),LSi_dic[s[0]+1]] = coeffs[fitorder-1] * self.__CONV         # in GPa unit 
                
                
                deltas = []
                deltasq = 0
                poly = np.poly1d(coeffs)
                
                for i in range(len(stress)):
                    delta = stress[i] - poly(strain[i])
                    deltas.append(delta)
                    deltasq += (delta)**2.0
                    
                self.__r[(emax,fitorder,LSi_dic[l])]=(np.sqrt(deltasq/len(strain)))
            
                self.__coeffs[(emax,fitorder,LSi_dic[l])] = (coeffs[fitorder-1])
                #print coeffs, fitorder
                if l==1: sigma1.append(coeffs[fitorder-1])
                elif l==2: sigma2.append(coeffs[fitorder-1])
                elif l==3: sigma3.append(coeffs[fitorder-1])
                elif l==4: sigma4.append(coeffs[fitorder-1])
                elif l==5: sigma5.append(coeffs[fitorder-1])
                elif l==6: sigma6.append(coeffs[fitorder-1])
                
                etacalc.append(emax)
                
                if (abs(strain[0]+emax) < 1.e-7):
                    strain.pop(0)
                    stress.pop(0)
                if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                    strain.pop()
                    stress.pop()
                
        for i in range(len(sigma1)): 
            
            sigma[str(round(etacalc[i],5))] = [sigma1[i], sigma2[i], sigma3[i], sigma4[i], sigma5[i], sigma6[i]] 
            print 'coeff ',sigma[str(round(etacalc[i],5))]
            
        self.__sigma = sigma    
        
                
    def get_2nd(self):
        return self.__Cij2nd
    
    def get_sigma(self):
        return self.__sigma
    
    def get_r(self):
        return self.__r
    
    def set_3rd(self, fitorder):
        """Evaluate 3rd order elastic constants from energy strain curves.
        
        Parameters
        ----------
        fitorder : integer 
            Order of polynomial energy-strain fit.
        """
        self.search_for_failed()
        
        self.__CONV = self.__vToGPa * math.factorial(3)*2.
        strain = copy(self.__strain)
        energy = copy(self.__energy)
        while (len(strain) > fitorder): 
            
            emax  = max(strain)
            emin  = min(strain)
            emax  = max(abs(emin),abs(emax))
            self.__coeffs= np.polyfit(strain, energy, fitorder)
            self.__Cij3rd  = self.__coeffs[fitorder-3]*self.__CONV/self.__V0 * 0.001 # in TPa units
            
            """  Calculate RMS:  """
            deltasq = 0
            deltas = []
            
            poly = np.poly1d(self.__coeffs)
            
            for i in range(len(energy)):
                delta = energy[i] - poly(strain[i])
                deltas.append(delta)
                deltasq += (delta)**2.0
                
            self.rms.append(np.sqrt(deltasq/len(strain)))
            
            if (abs(strain[0]+emax) < 1.e-7):
                strain.pop(0); energy.pop(0)
            if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                strain.pop()
                energy.pop()
                
        
    def get_3rd(self):
        return self.__Cij3rd
    
    
    def set_cvs(self, fitorder):
        """Evaluation of cross validation score.
        
        Parameters
        ----------
        fitorder : integer 
            Order of polynomial energy-strain fit.
        """
        self.search_for_failed()
        
        
        LSi_dic = {1:'LS1',2:'LS2',3:'LS3',4:'LS4',5:'LS5',6:'LS6'}
        for l in range(1,7):

            
            if l==1: s=(0,0)
            elif l==2: s=(1,1)  
            elif l==3: s=(2,2) 
            elif l==4: s=(1,2) 
            elif l==5: s=(0,2) 
            elif l==6: s=(0,1) 
            stress_ii = [st[s[0]][s[1]] for st in self.__stress]
            
        
            strain = copy(self.__strain)
            stress = copy(stress_ii)
            
            while (len(strain) > fitorder+1):
                emax = max(strain)
                emin = min(strain)
                emax = max(abs(emin),abs(emax))
        
                S = 0
                for k in range(len(strain)):
                    Y      = stress[k]
                    etatmp = []
                    enetmp = []
    
                    for h in range(len(strain)):
                        if (h==k): pass
                        else:
                            etatmp.append(strain[h])
                            enetmp.append(stress[h])
        
                    Yfit = np.polyval(np.polyfit(etatmp,enetmp, fitorder), strain[k])
                    S    = S + (Yfit-Y)**2
                
                self.__CV.append((np.sqrt(S/len(strain)),emax,fitorder))
                
        
                if (abs(strain[0]+emax) < 1.e-7):
                    strain.pop(0)
                    stress.pop(0)
                if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                    strain.pop()
                    stress.pop()
                
                
    def get_cvs(self):
        return self.__CV
    
    def plot_stress(self, etamax=0.05, fitorder=4):
        """Return matplotlib axis instance for energy-strain curve.  """
        if not mpl: raise "Problem with matplotib: Plotting not possible."
        self.search_for_failed()
        #f = plt.figure(figsize=(5,4), dpi=100)
        #ax = f.add_subplot(111)
        plt.plot(self.__strain, self.__energy, '*')
        
        poly = np.poly1d(self.__coeffs[(etamax,fitorder)])
        xp = np.linspace(min(self.__strain), max(self.__strain), 100)
        ax = plt.plot(xp, poly(xp))
        
        plt.xlabel('strain')
        plt.ylabel(r'energy    in eV')
        return ax
        
    def search_for_failed(self, mod='pass'):
        """Handle failed calculations"""
        i=0
        while i < len(self.__strain):
            
            if self.__stress[i] == 0:
                print "WARNING: Calculation failed on %s ...... trying to fit anyway."%str(self.__strain[i])
                self.__stress.pop(i)
                self.__strain.pop(i)
            i+=1
            
    stress = property( fget = get_stress        , fset = set_stress    )
    strain = property( fget = get_strain        , fset = set_strain    )
    V0 = property( fget = get_V0        , fset = set_V0    )
                
        
        