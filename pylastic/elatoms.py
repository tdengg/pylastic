'''
Created on May 23, 2014

@author: t.dengg
'''
import os
import copy
import cPickle as pickle
import numpy as np

from distort import Distort
from spacegroup import Sgroup
from vaspIO import POS

class ElAtoms(Distort, Sgroup, POS):
    '''
    ASE Atoms like object.
    
    
    Methods:
        set_cell
        set_natom
        set_scale
        set_species
        set_workdir
        set_poscar
        set_poscarnew
        distort
        poscarToAtoms
        atomsToPoscar
     
        get_cell
        get_natom
        get_scale
        get_species
        get_workdir
        get_poscar
        get_poscarnew
     
    
    Example:


    .. code-block:: python

        poscar = POS('POSCAR').read_pos()
        structures = Structures()
    
    
    Generate distortion:

    .. code-block:: python

        atom1 = ElAtoms()
        atom1.poscarToAtoms(poscar)
        atom1.distort(eta=0.01, strainType_index = 0)
        structures.append_structure(atom1)
    
    Generate another distortion:

    .. code-block:: python

        atom2 = ElAtoms()
        atom2.poscarToAtoms(poscar)
        atom2.distort(eta=0.01, strainType_index = 1)
        structures.append_structure(atom2)
    
        structures.write_structures()
    '''
    

    def __init__(self):
        '''
        Constructor
        '''
        Distort.__init__(self)
        
        self.__poscarnew = None
        self.__cell = None
        self.__gsenergy = None
        self.__executable = None
        self.__V0 = None
        
        
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
    
    def set_path(self, path):
        self.__path = path
        
    def get_path(self):
        return self.__path
    
    def distort(self, eta=0.0, strainType_index=0):
        self.sgn = self.__sgn
        
        self.strainType = self.get_strainList()[strainType_index]
        
        self.eta = eta
        self.set_defMatrix()
        def_matrix = self.defMatrix
        M_new = np.dot(self.__cell, def_matrix)
        self.__cell = M_new
        self.__poscarnew = copy.deepcopy(self.__poscar)
        self.__poscarnew['vlatt_1'] = self.__cell[0]
        self.__poscarnew['vlatt_2'] = self.__cell[1]
        self.__poscarnew['vlatt_3'] = self.__cell[2]
        
    def poscarToAtoms(self, poscar):
        self.__poscar = poscar
        self.__cell = np.array([self.__poscar['vlatt_1'],self.__poscar['vlatt_2'],self.__poscar['vlatt_3']])
        self.__natom = self.__poscar['natoms']
        self.__species = self.__poscar['vbasis']
        self.__scale = self.__poscar['scale']
        D = np.linalg.det(self.__cell)
        self.__V0 = abs(self.__scale**3*D)
        
        self.sgn = Sgroup(self.__poscar, self.__poscar['path']).sgn
        self.__sgn = self.sgn
    
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
    
    def set_poscarnew(self, poscar):
        self.__poscar = poscar
        
    def get_poscarnew(self):
        return self.__poscarnew
    
    def set_V0(self, V0):
        self.__V0 = V0
        
    def get_V0(self):
        return self.__V0
    
    def set_gsenergy(self, gsenergy):
        self.__gsenergy = gsenergy
        
    def get_gsenergy(self):
        return self.__gsenergy
    
    V0      = property( fget = get_V0       , fset = set_V0   )
    cell    = property( fget = get_cell         , fset = set_cell    )
    natom   = property( fget = get_natom        , fset = set_natom   )
    scale   = property( fget = get_scale        , fset = set_scale   )
    species = property( fget = get_species      , fset = set_species )
    path    = property( fget = get_path         , fset = set_path    )
    poscar  = property( fget = get_poscar       , fset = set_poscar  )
    poscarnew  = property( fget = get_poscarnew       , fset = set_poscarnew )
    gsenergy= property( fget = get_gsenergy     , fset = set_gsenergy)
    
    
class Structures(ElAtoms, Sgroup, POS):
    """Generate a series of distorted structures.
    
    Methods:
     set_fname
     write_structures
     append_structure
     get_structures
    
    """
    def __init__(self):
        ElAtoms.__init__(self)
        self.__structures = {}
        self.__fnames = []
        self.__workdir = '.'
        
    def set_fname(self, fnames):
        self.__fnames = fnames
        
    def set_workdir(self, workdir):
        self.__workdir = workdir
        
    def get_workdir(self):
        return self.__workdir
        
    def write_structures(self):
        dirnames = []
        os.chdir(self.__workdir)
        for atoms in self.__structures:
            try:
                if not atoms[0] in dirnames: os.mkdir(atoms[0]) 
                os.mkdir('%s/%s'%(atoms[0],atoms[1]))
            except:
                print 'dir exists'
            dirnames.append(atoms[0])
            fname = '%s/%s'%(atoms[0],atoms[1])
            self.__fnames.append(fname)
            os.system('cp KPOINTS %s/%s/KPOINTS'%(atoms[0],atoms[1]))
            os.system('cp INCAR %s/%s/INCAR'%(atoms[0],atoms[1]))
            os.system('cp POTCAR %s/%s/POTCAR'%(atoms[0],atoms[1]))
            POS().write_pos(self.__structures[atoms].poscarnew, fname+'/POSCAR')
            self.__structures[atoms].path = self.__workdir + fname
            obj = self.__structures[atoms]
            
            with open(fname+'/atoms.pkl', 'wb') as output:
                pickle.dump(obj, output, -1)
                
        with open('structures.pkl', 'wb') as output:
            pickle.dump(self.__structures, output, -1)
    
    def set_executable(self, executable):
        self.__executable = executable
        
    def get_executable(self):
        return self.__executable
    
    
    def calc_vasp(self):
        for atoms in self.__structures:
            
            fname = '%s/%s'%(atoms[0],atoms[1])
            self.__fnames.append(fname)
            os.chdir('%s/%s/'%(atoms[0],atoms[1]))
            os.system(self.__executable)
            os.chdir('../../')
        
    def append_structure(self, atoms):
        self.__structures[(atoms.strainType,round(atoms.eta,3))]= atoms
        
    def get_atomsByStraintype(self, strainType):
        atomslist = []
        for dic in sorted(self.__structures):
            if dic[0] == strainType: 
                atomslist.append(self.__structures[dic])
        return atomslist
        
        
    def get_structures(self):
        return self.__structures
    
    
    executable    = property( fget = get_executable        , fset = set_executable)
    workdir = property( fget = get_workdir        , fset = set_workdir)
    
    