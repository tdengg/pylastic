'''
'''
import os
import copy
import time

from subprocess import Popen, PIPE
import cPickle as pickle
import numpy as np

from pylastic.distort import Distort
from pylastic.spacegroup import Sgroup
from pylastic.io.vasp import POS
from pylastic.prettyPrint import PrettyMatrix
from pylastic.status import Check

class ElAtoms(Distort, Sgroup, POS, PrettyMatrix, Check):
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
        super(Check, self).__init__()
        
        self.__poscarnew = None
        self.__cell = None
        self.__gsenergy = None
        self.__executable = None
        self.__V0 = None
        self.__verbose = True
        
        self.__mthd = 'Energy'
        
        
    def set_cell(self, cell):
        """Set lattice vectors of crystal cell."""
        if isinstance(cell, list): self.__cell = cell
        else: print 'Invalide Type of lattice vector: %s. Must be list instead!'%(type(cell))
        
        
    def get_cell(self):
        return self.__cell
    
    
    def set_natom(self, natom):
        """Set number of atoms in supercell.
        
        Parameters
        ----------
        natom : integer
            Supercell size.
        """
        if isinstance(natom, int): self.__natom = natom
        else: print 'Number of atoms is invalid type: %s. Must be integer instead!'%(type(natom))
        
    def get_natom(self):
        return self.__natom
    
    
    def set_scale(self, scale):
        """Set scaling factor for cell vectors
        
        Parameters
        ----------
        scale : float
            Unit cell scaling.
        """
        if isinstance(scale, float) or isinstance(scale, list): self.__scale = scale
        else: print 'Scale is of invalid type: %s. Must be float (or list of float values) instead!'%(type(scale))
        
    def get_scale(self):
        return self.__scale
    
    
    def set_species(self,species):
        """Set basis vectors.
        
        Parameters
        ----------
        species : list
            Basis vectors
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
        self.mthd = self.__mthd
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
    
    
        
    ### VASP
    ### --->    
    def poscarToAtoms(self, poscar):
        """Set properties of poscar input to atoms object.
        
        Parameters
        ----------
        poscar : dictionary
            POSCAR file converted to dictionary.
        """
        if self.__verbose: print 'Converting dictionary to atoms object ....'
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
    ### <---
    
    
    
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
    
    
    def set_method(self, mthd):
        if mthd in ['Energy','Stress']: self.__mthd = mthd
        else: print "Wrong value for method: Please choose either 'Energy' or 'Stress'!"
    
    def get_method(self):
        return self.__mthd
    
    
    V0      = property( fget = get_V0       , fset = set_V0   )
    cell    = property( fget = get_cell         , fset = set_cell    )
    natom   = property( fget = get_natom        , fset = set_natom   )
    scale   = property( fget = get_scale        , fset = set_scale   )
    species = property( fget = get_species      , fset = set_species )
    path    = property( fget = get_path         , fset = set_path    )
    poscar  = property( fget = get_poscar       , fset = set_poscar  )
    poscarnew  = property( fget = get_poscarnew       , fset = set_poscarnew )
    gsenergy= property( fget = get_gsenergy     , fset = set_gsenergy)
    mthd = property( fget = get_method       , fset = set_method)
    
    
    
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
        if os.path.isdir(workdir): self.__workdir = workdir
        else: "Directory '%s' does not exist!"%workdir
        
    def get_workdir(self):
        return self.__workdir
        
    def write_structures(self, strkt, overwrite=True):
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
            
        Parameters
        ----------
        strkt : object
            Structures object to be written to structures.pkl
        
        Keyword arguments
        -----------------
        overwrite : boolean
            Specify if existing files should be overwritten (default=True)
        """
        dirnames = []
        os.chdir(self.__workdir)
        for atoms in self.__structures:
            try:
                if not atoms[0] in dirnames: os.mkdir(atoms[0]) 
                os.mkdir('%s/%s'%(atoms[0],atoms[1]))
            except:
                print "Directory '%s/%s' already exists."%(atoms[0],atoms[1])
            dirnames.append(atoms[0])
            fname = '%s/%s'%(atoms[0],atoms[1])
            self.__fnames.append(fname)
            
            if not os.path.isfile("%s/%s/KPOINTS"%(atoms[0],atoms[1])) or overwrite: os.system('cp KPOINTS %s/%s/KPOINTS'%(atoms[0],atoms[1]))
            else: print "%s/%s/KPOINTS already existing: overwrite = False"%(atoms[0],atoms[1])
            if not os.path.isfile("%s/%s/INCAR"%(atoms[0],atoms[1]))   or overwrite: os.system('cp INCAR %s/%s/INCAR'%(atoms[0],atoms[1]))
            else: print "%s/%s/INCAR   already existing: overwrite = False"%(atoms[0],atoms[1])
            if not os.path.isfile("%s/%s/POTCAR"%(atoms[0],atoms[1]))  or overwrite: os.system('cp POTCAR %s/%s/POTCAR'%(atoms[0],atoms[1]))
            else: print "%s/%s/POTCAR  already existing: overwrite = False"%(atoms[0],atoms[1])
            if not os.path.isfile(fname+'/POSCAR')                     or overwrite: POS().write_pos(self.__structures[atoms].poscarnew, fname+'/POSCAR')
            else: print "%s/%s/POSCAR  already existing: overwrite = False"%(atoms[0],atoms[1])
            self.__structures[atoms].path = self.__workdir + fname
            obj = self.__structures[atoms]
            
            with open(fname+'/atoms.pkl', 'wb') as output:
                pickle.dump(obj, output, -1)
                
        with open('structures.pkl', 'wb') as output:
            pickle.dump(strkt, output, -1)
        
    
    def set_executable(self, executable):
        """Set path to vasp executable. 
        
        Parameters
        ----------
        executable : string 
            Location of vasp executable (e.g. /home/user/bin/vasp)
        """
        
        if os.path.isfile(executable): self.__executable = executable
        else: "Wrong path to executable '%s'"%executable
        
    def get_executable(self):
        return self.__executable
    
    
    def calc_vasp(self,lock=None, overwrite = True):
        """Perform local vasp calculations.
        
        Keyword arguments
        -----------------
        lock : object 
            thread locker object (default=None)
        overwrite : boolean
            Specify if existing files should be overwritten (default=True)
        """
        
        for atoms in self.__structures:
            
            fname = '%s/%s'%(atoms[0],atoms[1])
            self.__fnames.append(fname)
            os.chdir('%s/%s/'%(atoms[0],atoms[1]))
            string = "# Starting vasp calculation for %s/%s ......... #"%(atoms[0],atoms[1])
            n = len(string)
            
            print "\n"+ n*'#'
            print string
            print n*'#' + "\n"
            if not lock==None: lock.acquire()
            sp = Popen([self.__executable])
            sp.communicate()
            time.sleep(0.01)
            if not lock==None: lock.release()
            
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
    
    
    def status(self):
        state = Check()
        state.workdir = self.__workdir
        self.__status = state.check_calc()
        return self.__status
        
    executable    = property( fget = get_executable        , fset = set_executable)
    workdir = property( fget = get_workdir        , fset = set_workdir)
    
    
    
    