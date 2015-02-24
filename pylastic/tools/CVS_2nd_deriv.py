import numpy as np
from copy import copy

class CVS(object):
    def __init__(self, xdata, ydata):
        self.__strain = list(xdata)
        self.__energy = list(ydata)
    
    def calc_cvs(self, fitorder):
        
        
        
        strain = copy(self.__strain)
        energy = copy(self.__energy)
        
        self.__cvs = []
        
        while (len(strain) > fitorder+1):
            emax = max(strain)
            emin = min(strain)
            emax = max(abs(emin),abs(emax))
    
            S = 0
            for k in range(len(strain)):
                #Y      = energy[k]
                Y = np.polyfit(strain,energy, fitorder)[fitorder-2]*2.
                etatmp = []
                enetmp = []

                for l in range(len(strain)):
                    if (l==k): pass
                    else:
                        etatmp.append(strain[l])
                        enetmp.append(energy[l])
                
                Yfit = np.polyfit(etatmp,enetmp, fitorder)[fitorder-2]*2.
                #Yfit = np.polyval(np.polyfit(etatmp,enetmp, fitorder), strain[k])
                S    = S + (Yfit-Y)**2
            
            self.__cvs.append((np.sqrt(S/len(strain)),emax,fitorder))
            
    
            if (abs(strain[0]+emax) < 1.e-7):
                strain.pop(0)
                energy.pop(0)
            if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                strain.pop()
                energy.pop()
                
        
        
        return self.__cvs
    
    