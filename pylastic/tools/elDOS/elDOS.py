from scipy.integrate import simps
import matplotlib.pyplot as plt
import numpy as np

class DOS(object):
    def __init__(self, D=None):
        self.__D = D
    
    def read_DOSCAR(self, fname):
        f=open(fname)
        lines = f.readlines()
        f.close()
        
        energy = [float(l.split()[0]) for l in lines[6:]]
        D = [float(l.split()[1]) for l in lines[6:]]
        sumD = [float(l.split()[2]) for l in lines[6:]]
        print len(energy),len(D)
        self.__D = D
        self.__energy = energy
    
    def plot_DOS(self):
        plt.plot(self.__energy,self.__D)
        
if __name__ == "__main__":
    dos = DOS()
    dos.read_DOSCAR("DOSCAR")
    dos.plot_DOS()
    plt.show()