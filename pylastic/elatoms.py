'''
Created on May 23, 2014

@author: t.dengg
'''
import os
import numpy as np
import pickle
from distort import Distort
from spacegroup import Sgroup
import numpy as np
from vaspIO import POS

class ElAtoms(Distort, Sgroup, POS):
    '''
    classdocs:
    Generate ASE Atoms like object.
    '''
    

    def __init__(self):
        '''
        Constructor
        '''
        Distort.__init__(self)
        
        
        self.__cell = None
        
    def set_cell(self, cell):
        self.__cell = cell
        
    def get_cell(self):
        return self.__cell
    
    def set_natom(self, natom):
        self.__natom = natom
        
    def get_natom(self):
        return self.__natom
    
    def set_scale(self, scale):
        self.__scale = scale
        
    def get_scale(self):
        return self.__scale
    
    def set_species(self,species):
        self.__species = species
    
    def get_species(self):
        return self.__species
    
    def set_workdir(self, wdir):
        self.__workdir = wdir
        
    def get_workdir(self):
        return self.__workdir
    
    def distort(self, eta=None, strainType_index=None):
        self.sgn = self.__sgn
        if strainType_index:
            self.strainType = self.get_strainList()[strainType_index]
        else:
            self.strainType = next(self.get_strainList_iter())
            
        if eta:
            self.eta = eta
        self.set_defMatrix()
        def_matrix = self.defMatrix
        M_new = np.dot(self.__cell, def_matrix)
        self.__cell = M_new
    
    def poscarToAtoms(self, poscar):
        self.__poscar = poscar
        self.__cell = np.array([self.__poscar['vlatt_1'],self.__poscar['vlatt_2'],self.__poscar['vlatt_3']])
        self.__natom = self.__poscar['natoms']
        self.__species = self.__poscar['vbasis']
        self.__scale = self.__poscar['scale']
        
        self.__sgn = Sgroup(self.__poscar, self.__poscar['path']).sgn
        
    
    def atomsToPoscar(self):
        self.__poscar = {}
        self.__poscar['vlatt_1'] = self.__cell[0]
        self.__poscar['vlatt_2'] = self.__cell[1]
        self.__poscar['vlatt_3'] = self.__cell[2]
        self.__poscar['natoms'] = self.__natom
        self.__poscar['vbasis'] = self.__species
        self.__poscar['scale'] = self.__scale
        
        return self.__poscar
        
    
    def set_poscar(self, poscar):
        self.__poscar = poscar
        
    def get_poscar(self):
        return self.__poscar
    
    cell    = property( fget = get_cell         , fset = set_cell   )
    natom   = property( fget = get_natom        , fset = set_natom  )
    scale   = property( fget = get_scale        , fset = set_scale  )
    species = property( fget = get_species      , fset = set_species)
    workdir = property( fget = get_workdir      , fset = set_workdir)
    poscar  = property( fget = get_poscar       , fset = set_poscar )
    
    
class Structures(ElAtoms, Sgroup, POS):
    """Generate a series of distorted atoms.
    """
    def __init__(self):
        ElAtoms.__init__(self)
        
        self.__structures = []
        self.__fname = os.getcwd()+'/POSCAR'
        
    def set_fname(self, fname):
        self.__fname = fname
        
    def write_structures(self):
        for atoms in self.__structures:
            POS.write_pos(atoms.poscar, self.__fname)
            
    def append_structure(self, atoms):
        self.__structures.append(atoms)
        
    def get_structures(self):
        return self.__structures
    