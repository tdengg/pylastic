import numpy as np

class CALC(object):
    def __init__(self, verbous=True):
        self.__verbous = verbous
        self.__code = 'vasp'
    def read_POS(self, fname):
        if self.__code == 'vasp':
            from pylastic.io.vasp import POS
            poscar = POS().read_pos()
        self.__poscar = poscar
        print self.__poscar
    
    def set_POS(self, poscar):
        self.__poscar = poscar
    
    def get_POS(self):
        return self.__poscar
    
    def set_a(self,a):
        self.__a = a
    def get_a(self):
        return self.__a
    
    def set_b(self,b):
        self.__b = b
    def get_b(self):
        return self.__b
    
    def set_c(self,c):
        self.__c = c
    def get_c(self):
        return self.__c
        
    def volume(self):
        scale = self.__poscar['scale']
        if scale<0.:
            volumemod=True
            vol = -scale
            scale=1.
        else:
            volumemod=False
        self.__scale = scale
        self.__cell = np.array([self.__poscar['vlatt_1'],self.__poscar['vlatt_2'],self.__poscar['vlatt_3']])*self.__scale
        self.__a = np.sqrt(self.__poscar['vlatt_1'][0]**2. + self.__poscar['vlatt_1'][1]**2. + self.__poscar['vlatt_1'][2]**2.)*self.__scale
        self.__b = np.sqrt(self.__poscar['vlatt_2'][0]**2. + self.__poscar['vlatt_2'][1]**2. + self.__poscar['vlatt_2'][2]**2.)*self.__scale
        self.__c = np.sqrt(self.__poscar['vlatt_3'][0]**2. + self.__poscar['vlatt_3'][1]**2. + self.__poscar['vlatt_3'][2]**2.)*self.__scale
        
        natoms = 0
        print self.__cell, self.__scale
        for nion in self.__poscar['natoms']:
            natoms += nion
        if volumemod:
            V = vol
        else:
            V = np.linalg.det(self.__cell)/natoms
        
        return V
    
    poscar = property( fget = get_POS       , fset = set_POS)
    a = property( fget = get_a       , fset = set_a)
    b = property( fget = get_b       , fset = set_b)
    c = property( fget = get_c       , fset = set_c)