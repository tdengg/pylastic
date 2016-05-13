import numpy as np
import os

class Generate_XSF(object):
    def __init__(self, fname_outcar):
        self.__fname_outcar= fname_outcar
        
    
    def read_eigenvect_OUTCAR(self):
        f = open(self.__fname_outcar)
        lines=f.readlines()
        for line in lines:
            if 'Eigenvectors' in line.split():
                
            else:
                continue