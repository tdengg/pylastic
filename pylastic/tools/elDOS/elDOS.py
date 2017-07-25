from scipy.integrate import simps
import matplotlib.pyplot as plt
import numpy as np
import re

class DOS(object):
    def __init__(self, D=None):
        self.__k_B = 8.6173325*10**(-5)#eV/K
        self.__D = D
    
    def read_DOSCAR(self, path):
        f=open(path+'/DOSCAR')
        lines = f.readlines()
        f.close()
        
        energy = [float(l.split()[0]) for l in lines[7:]]
        D = [float(l.split()[1]) for l in lines[7:]]
        sumD = [float(l.split()[2]) for l in lines[7:]]
        
        self.__D = D
        self.__energy = energy
    
    def plot_DOS(self):
        plt.plot(self.__energy,self.__D)
        
    def free_energy(self, T):
        
        fd = lambda engy: 1./( 1.+np.exp((engy-self.__efermi)/(self.__k_B*T)) )
        
        integrand = []
        for e,d in zip(self.__energy, self.__D):
            
            if (e-self.__efermi)/(self.__k_B*T) <= 50000.:
                integrand.append(-self.__k_B*T*d*fd(e)*np.log(fd(e))*((1.-fd(e))*(1.-np.log(fd(e)))))
            else:
                integrand.append(0.)
        #plt.plot(self.__energy,[fd(engy)for engy in self.__energy])
        
        return simps(integrand)
    
    def get_eFermi(self,path):
        with open(path+"/OUTCAR") as origin_file:
            for line in origin_file:
                if 'E-fermi' in line:
            
                    self.__efermi = float(line.split()[2])
        return self.__efermi
    
if __name__ == "__main__":
    dos = DOS()
    dos.read_DOSCAR("DOSCAR")
    dos.plot_DOS()
    
    dos.free_energy(1000.)
    
    plt.show()