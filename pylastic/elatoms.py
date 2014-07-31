'''
'''
import os
import copy
import cPickle as pickle
import numpy as np

from distort import Distort
from spacegroup import Sgroup
from vaspIO import POS
from prettyPrint import PrettyMatrix

class ElAtoms(Distort, Sgroup, POS, PrettyMatrix):
    '''
    ASE Atoms like object.
     
    **Example:**

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
        super(ElAtoms, self).__init__()
        PrettyMatrix.__init__(self)
        
        self.__poscarnew = None
        self.__cell = None
        self.__gsenergy = None
        self.__executable = None
        self.__V0 = None
        
        
    def set_cell(self, cell):
        """Set lattice vectors of crystal cell."""
        self.__cell = cell
        
        
    def get_cell(self):
        return self.__cell
    
    def set_natom(self, natom):
        """Set number of atoms in supercell.
        
        Parameters
        ----------
        natom : integer
            Supercell size.
        """
        self.__natom = natom
        
    def get_natom(self):
        return self.__natom
    
    def set_scale(self, scale):
        """Set scaling factor for cell vectors
        
        Parameters
        ----------
        scale : float
            Unit cell scaling.
        """
        self.__scale = scale
        
    def get_scale(self):
        return self.__scale
    
    def set_species(self,species):
        """Set atoms species.
        
        Parameters
        ----------
        species : string
            Atom species
        """
        self.__species = species
    
    def get_species(self):
        return self.__species
    
    def set_path(self, path):
        """Set path to the calculations sub-directory.
        
        Parameters
        ----------
        path : string
            Path to calculation subdir.
        """
        self.__path = path
        
    def get_path(self):
        return self.__path
    
    def distort(self, eta=0.0, strainType_index=0):
        """Distort structure and 
        
        Parameters
        ----------
        eta : float
            lagrangian strain value
            
        strainType_index : integer
            List-index of strainType in strainList. 
        """
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
        """Set properties of poscar input to atoms object.
        
        Parameters
        ----------
        poscar : dictionary
            POSCAR file converted to dictionary.
        """
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
        """Create poscar out of atoms object."""
        self.__poscar = {}
        self.__poscar['vlatt_1'] = self.__cell[0]
        self.__poscar['vlatt_2'] = self.__cell[1]
        self.__poscar['vlatt_3'] = self.__cell[2]
        self.__poscar['natoms'] = self.__natom
        self.__poscar['vbasis'] = self.__species
        self.__poscar['scale'] = self.__scale
        
        return self.__poscar
    
    def set_poscar(self, poscar):
        """"""
        self.__poscar = poscar
        
    def get_poscar(self):
        return self.__poscar
    
    def set_poscarnew(self, poscar):
        """"""
        self.__poscar = poscar
        
    def get_poscarnew(self):
        return self.__poscarnew
    
    def set_V0(self, V0):
        """Set equilibrium volume of supercell.
        
        Parameters
        ----------
        V0 : float
            Equilibrium volume of supercell.
        """
        self.__V0 = V0
        
    def get_V0(self):
        return self.__V0
    
    def set_gsenergy(self, gsenergy):
        """Give atoms object a groundstate energy.
        
        Parameters
        ----------
        gsenergy : float
            groundstate energy in eV
        """
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
    
    """
    def __init__(self):
        super(Structures, self).__init__()
        self.__structures = {}
        self.__fnames = []
        self.__workdir = './'
        
    def set_fname(self, fnames):
        """Set absolute path to calculations sub-directory.
        
        Parameters
        ----------
        fnames : string
            Absolute pathname to calculation sub-directory.
        """
        self.__fnames = fnames
        
    def set_workdir(self, workdir):
        """Set working directory.
        
        Parameters
        ----------
        workdir : string
            Working directory.
        """
        self.__workdir = workdir
        
    def get_workdir(self):
        return self.__workdir
        
    def write_structures(self):
        """Generate a file structure and write all input files for vasp. 
        
        Filestructure:
        
        * root
            * straintype1
                * eta1
                * eta2
                * ...
            * straintype2
                * eta1
                * eta2
                * ...
            * ...
        """
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
        """Set path to vasp executable. 
        
        Parameters
        ----------
        executable : string 
            Location of vasp executable (e.g. /home/user/bin/vasp)
        ."""
        self.__executable = executable
        
    def get_executable(self):
        return self.__executable
    
    
    def calc_vasp(self):
        """Perform local vasp calculations."""
        for atoms in self.__structures:
            
            fname = '%s/%s'%(atoms[0],atoms[1])
            self.__fnames.append(fname)
            os.chdir('%s/%s/'%(atoms[0],atoms[1]))
            os.system(self.__executable)
            os.chdir('../../')
        
    def append_structure(self, atoms):
        """Append structure to Structures object to create a set of distorted structures.
        
        Parameters
        ----------
        atoms : object 
            atoms object created from ElAtoms class.
        """
        self.__structures[(atoms.strainType,round(atoms.eta,3))]= atoms
        
    def get_atomsByStraintype(self, strainType):
        """Returns a sorted list of structures with given strain-type.
        
        Parameters
        ----------
        strainType : string 
            Strain type.
        """
        atomslist = []
        for dic in sorted(self.__structures):
            if dic[0] == strainType: 
                atomslist.append(self.__structures[dic])
        return atomslist
        
        
    def get_structures(self):
        """Return dictionary with all structures."""
        return self.__structures
    
    
    executable    = property( fget = get_executable        , fset = set_executable)
    workdir = property( fget = get_workdir        , fset = set_workdir)
    
    