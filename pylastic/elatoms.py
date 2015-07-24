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

from pylastic.prettyPrint import PrettyMatrix
from pylastic.status import Check

class ElAtoms(Distort, Sgroup, PrettyMatrix, Check):
    '''
    Object similar to the ``atoms`` class in *ASE*, containing information related to the structure.
     
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
    

    def __init__(self, cod='vasp'):
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
        self.__code = cod
        self.__path = os.getcwd()
        
        
        
    def set_cell(self, cell):
        """Set lattice vectors of crystal cell.
        
        cell : list of 3 lists
            Lattice vectors in Cartesian coordinates.
        """
        
        if isinstance(cell, list): self.__cell = cell
        else: print 'Invalide Type of lattice vector: %s. Must be list instead!'%(type(cell))
        
        
    def get_cell(self):
        return self.__cell
    
    
    def set_natom(self, natom):
        """Set number of atoms in crystal cell.
        
        
        
        natom : integer
            Number of atoms in crystal cell.
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
        if isinstance(scale, float) or isinstance(scale, list): 
            self.__scale = scale
            self.__V = np.linalg.det(self.__cell*scale)
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
        """Distort the structure.
        
        Parameters
        ----------
        eta : float
            lagrangian strain value
            
        strainType_index : integer
            List-index of strainType in strainList. 
        """
        #self.mthd = self.__mthd
        self.sgn = self.__sgn
        if not strainType_index==-1:
            self.strainType = self.get_strainList()[strainType_index]
        else:
            self.strainType = 'vol'
        self.eta = eta
        self.set_defMatrix(self.__code)
        def_matrix = self.defMatrix
        M_new = np.dot(self.__cell, def_matrix)
        self.__cell = M_new
        self.__V = np.linalg.det(self.__cell*self.__scale)
        self.__poscarnew = copy.deepcopy(self.__poscar)
        self.__poscarnew['vlatt_1'] = self.__cell[0]
        self.__poscarnew['vlatt_2'] = self.__cell[1]
        self.__poscarnew['vlatt_3'] = self.__cell[2]
        
    
    def deform_volume(self, eta=0.0):
        """Hydrostatic deformation"""    
        def_matrix = np.array([[1,0,0],[0,1,0],[0,0,1]])+np.array([[eta,0,0],[0,eta,0],[0,0,eta]])
        M_new = np.dot(self.__cell, def_matrix)
        self.eta=eta
        self.__cell = M_new
        self.__poscarnew = copy.deepcopy(self.__poscar)
        self.__poscarnew['vlatt_1'] = self.__cell[0]
        self.__poscarnew['vlatt_2'] = self.__cell[1]
        self.__poscarnew['vlatt_3'] = self.__cell[2]
        #print np.linalg.det(self.__cell*self.__scale)
        self.__V = np.linalg.det(self.__cell*self.__scale)
        self.set_strainType('vol')
        
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
        if self.__code=='vasp': self.sgn = Sgroup(self.__poscar, self.__poscar['path']).sgn
        else: self.sgn = self.__poscar['sgn']
        self.__sgn = self.sgn
    
    def atomsToPoscar(self):
        """Create poscar out of atoms object."""
        poscar = {}
        print 'creating poscar: scale=%s'%self.__scale
        poscar['vlatt_1'] = self.__cell[0]
        poscar['vlatt_2'] = self.__cell[1]
        poscar['vlatt_3'] = self.__cell[2]
        poscar['natoms'] = self.__natom
        poscar['vbasis'] = self.__species
        poscar['scale'] = self.__scale
        ### Not defined
        poscar['selective'] = False
        poscar['csystem'] = 'direct'
        return poscar
    
    def set_poscar(self, poscar):
        """Assign the dictionary poscar."""
        self.__poscar = poscar
        
    def get_poscar(self):
        """Assign the dictionary poscar."""
        return self.__poscar
    
    def set_poscarnew(self, poscar):
        """"""
        self.__poscarnew = poscar
        
    def get_poscarnew(self):
        return self.__poscarnew
    ### <---
    
    
    
    def set_V0(self, V0):
        """Set equilibrium volume of crystal cell.
        
        Parameters
        ----------
        V0 : float
            Equilibrium volume of crystal cell.
        """
        self.__V0 = V0
        
    def get_V0(self):
        return self.__V0
    
    def set_V(self, V):
        """Set volume of distorted cell.
        
        Parameters
        ----------
        V : float
            Volume of distorted crystal cell.
        """
        self.__V = V
        
    def get_V(self):
        return self.__V
    
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
    
    def set_phenergy(self, phenergy):
        """Give atoms object a phonon free energy.
        
        Parameters
        ----------
        phenergy : float
            groundstate energy in eV
        """
        self.__phenergy = phenergy
        
    def get_phenergy(self):
        return self.__phenergy
    
#    def set_method(self, mthd):
#        if mthd in ['Energy','Stress']: self.__mthd = mthd
#        else: print "Wrong value for method: Please choose either 'Energy' or 'Stress'!"
#    
#    def get_method(self):
#        return self.__mthd
    
    def set_code(self, code):
        if code in ['vasp','exciting','espresso','wien']:
            self.__code = code
        else:
            print "Unknown code '%s'. Please choose either espresso, exciting, wien or vasp"%code
            
    def get_code(self):
        return self.__code
    
    def set_T(self,T):
        self.__T = T
    def get_T(self):
        return self.__T

    
    
    V0      = property( fget = get_V0       , fset = set_V0   )
    V       = property( fget = get_V      , fset = set_V   )
    cell    = property( fget = get_cell         , fset = set_cell    )
    natom   = property( fget = get_natom        , fset = set_natom   )
    scale   = property( fget = get_scale        , fset = set_scale   )
    species = property( fget = get_species      , fset = set_species )
    path    = property( fget = get_path         , fset = set_path    )
    poscar  = property( fget = get_poscar       , fset = set_poscar  )
    poscarnew  = property( fget = get_poscarnew       , fset = set_poscarnew )
    gsenergy= property( fget = get_gsenergy     , fset = set_gsenergy)
    phenergy= property( fget = get_phenergy     , fset = set_phenergy)
    T = property( fget = get_T       , fset = set_T   )
    
#    mthd = property( fget = get_method       , fset = set_method)
    code = property( fget = get_code       , fset = set_code)
    
    
class Structures(ElAtoms, Sgroup):
    """Generate a series of distorted structures.
    
    """
    def __init__(self, cod = 'vasp', thermo=False):
        super(Structures, self).__init__()
        self.__structures = {}
        self.__fnames = []
        self.__workdir = os.getcwd()+'/'
        self.__code = cod
        self.__thermo = thermo
        
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
        
    def write_structures(self, struct, overwrite=True):
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
        struct : object
            ``Structures`` object to be written to structures.pkl
        
        Keyword arguments
        -----------------
        overwrite : boolean
            Specify if existing files should be overwritten (default=True)
        """
        
        if self.__thermo:
            dirnames = []
            pardirnames = []
            os.chdir(self.__workdir)
            for atoms in self.__structures:
                if self.__structures[atoms].poscarnew==None: self.__structures[atoms].poscarnew=self.atomsToPoscar()
                self.__path = '%s/%s/eta%s'%(atoms[2],atoms[0],atoms[1])
                self.__structures[atoms].path = '%s/%s/eta%s'%(atoms[2],atoms[0],atoms[1])
                
                try:
                    
                    tempp=str(atoms[2])+'/'+atoms[0]
                    if not str(atoms[2]) in pardirnames: os.mkdir(str(atoms[2]))
                    if not tempp in dirnames: os.mkdir(tempp) 
                    os.mkdir('%s'%(self.__path))
                except:
                    print "Directory '%s' already exists."%(atoms[0])
                dirnames.append(str(atoms[2])+'/'+atoms[0])
                pardirnames.append(str(atoms[2]))
                
                self.__fnames.append(self.__path)
                
                ##### VASP #####
                if self.__code=='vasp':
                    
                    from pylastic.io.vasp import POS
                    if not os.path.isfile("%s/KPOINTS"%(self.__path)) or overwrite: os.system('cp KPOINTS %s/KPOINTS'%(self.__path))
                    else: print "%s/KPOINTS already existing: overwrite = False"%(self.__path)
                    if not os.path.isfile("%s/INCAR"%(self.__path))   or overwrite: os.system('cp INCAR %s/INCAR'%(self.__path))
                    else: print "%s/INCAR   already existing: overwrite = False"%(self.__path)
                    if not os.path.isfile("%s/POTCAR"%(self.__path))  or overwrite: os.system('cp POTCAR %s/POTCAR'%(self.__path))
                    else: print "%s/POTCAR  already existing: overwrite = False"%(self.__path)
                    if not os.path.isfile(self.__path+'/POSCAR')                     or overwrite: POS().write_pos(self.__structures[atoms].poscarnew, self.__path+'/POSCAR')
                    else: print "%s/POSCAR  already existing: overwrite = False"%(self.__path)
                ################
                
                
                ### ESPRESSO ###
                if self.__code=='espresso':
                    from pylastic.io.espresso import POS
                    if not os.path.isfile(self.__path+'/ElaStic_PW.in') or overwrite: POS('ElaStic_PW.in').write_in(self.__structures[atoms].poscarnew, self.__path+'/ElaStic_PW.in')
                    else: print "%s/ElaStic_PW.in  already existing: overwrite = False"%(self.__path)
                ################
                
                
                ### EXCITING ###
                if self.__code=='exciting':
                    from pylastic.io.exciting import POS
                    if not os.path.isfile(self.__path+'/input.xml') or overwrite: POS('input.xml').write_in(self.__structures[atoms].poscarnew, self.__path+'/input.xml')
                    else: print "%s/input.xml already existing: overwrite = False"%(self.__path)
                ################
                
                ### WIEN2K ###
                if self.__code=='wien':
                    
                    from pylastic.io.wien import POS
                    if not os.path.isfile(self.__path+'/distorted_P.struct') or overwrite: POS(self.__structures[atoms].poscarnew['path']).write_in(self.__structures[atoms].poscarnew, self.__path)
                    else: print "%s/ already existing: overwrite = False"%(self.__path)
                
                ################
                self.__structures[atoms].path = self.__workdir + self.__path
                obj = self.__structures[atoms]
                
                with open(self.__path+'/atoms.pkl', 'wb') as output:
                    pickle.dump(obj, output, -1)
                    
            with open('structures.pkl', 'wb') as output:
                pickle.dump(struct, output, -1)
            
        
        else:
            dirnames = []
            os.chdir(self.__workdir)
            print self.__structures
            for atoms in self.__structures:
                if self.__structures[atoms].poscarnew==None: 
                    self.__structures[atoms].poscarnew=self.__structures[atoms].atomsToPoscar()
                    print self.__structures[atoms].poscarnew
                    print self.__structures[atoms].scale
                    print self.__structures[atoms].cell
                if atoms[0] != None:
                    self.__path = '%s/eta%s'%(atoms[0],atoms[1])
                    self.__structures[atoms].path = '%s/eta%s'%(atoms[0],atoms[1])
                else:
                    self.__path = 'scale_%s'%(self.__structures[atoms].scale)
                    self.__structures[atoms].path = 'scale_%s'%(self.__structures[atoms].scale)
                try:
                    if not atoms[0] in dirnames and atoms[0] != None: os.mkdir(atoms[0]) 
                    os.mkdir('%s'%(self.__path))
                except:
                    print "Directory '%s' already exists."%(atoms[0])
                dirnames.append(atoms[0])
                
                self.__fnames.append(self.__path)
                
                ##### VASP #####
                if self.__code=='vasp':
                    
                    from pylastic.io.vasp import POS
                    if not os.path.isfile("%s/KPOINTS"%(self.__path)) or overwrite: os.system('cp KPOINTS %s/KPOINTS'%(self.__path))
                    else: print "%s/KPOINTS already existing: overwrite = False"%(self.__path)
                    if not os.path.isfile("%s/INCAR"%(self.__path))   or overwrite: os.system('cp INCAR %s/INCAR'%(self.__path))
                    else: print "%s/INCAR   already existing: overwrite = False"%(self.__path)
                    if not os.path.isfile("%s/POTCAR"%(self.__path))  or overwrite: os.system('cp POTCAR %s/POTCAR'%(self.__path))
                    else: print "%s/POTCAR  already existing: overwrite = False"%(self.__path)
                    if not os.path.isfile(self.__path+'/POSCAR')                     or overwrite: POS().write_pos(self.__structures[atoms].poscarnew, self.__path+'/POSCAR')
                    else: print "%s/POSCAR  already existing: overwrite = False"%(self.__path)
                ################
                
                
                ### ESPRESSO ###
                if self.__code=='espresso':
                    from pylastic.io.espresso import POS
                    if not os.path.isfile(self.__path+'/ElaStic_PW.in') or overwrite: POS('ElaStic_PW.in').write_in(self.__structures[atoms].poscarnew, self.__path+'/ElaStic_PW.in')
                    else: print "%s/ElaStic_PW.in  already existing: overwrite = False"%(self.__path)
                ################
                
                
                ### EXCITING ###
                if self.__code=='exciting':
                    from pylastic.io.exciting import POS
                    if not os.path.isfile(self.__path+'/input.xml') or overwrite: POS('input.xml').write_in(self.__structures[atoms].poscarnew, self.__path+'/input.xml')
                    else: print "%s/input.xml already existing: overwrite = False"%(self.__path)
                ################
                
                ### WIEN2K ###
                if self.__code=='wien':
                    
                    from pylastic.io.wien import POS
                    if not os.path.isfile(self.__path+'/distorted_P.struct') or overwrite: POS(self.__structures[atoms].poscarnew['path']).write_in(self.__structures[atoms].poscarnew, self.__path)
                    else: print "%s/ already existing: overwrite = False"%(self.__path)
                
                ################
                self.__structures[atoms].path = self.__workdir + self.__path
                obj = self.__structures[atoms]
                
                with open(self.__path+'/atoms.pkl', 'wb') as output:
                    pickle.dump(obj, output, -1)
                    
            with open('structures.pkl', 'wb') as output:
                pickle.dump(struct, output, -1)
            
        
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
            
            fname = '%s'%(self.__path)
            self.__fnames.append(fname)
            os.chdir('%s/'%(self.__structures[atoms].path))
            print os.getcwd()
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
            print 'WORKDIR: ', self.workdir
            os.chdir(self.__workdir)
        
        
    def append_structure(self, atoms):
        """Append structure to ``Structures`` object to create a set of distorted structures.
        
        Parameters
        ----------
        atoms : object 
            atoms object created from ElAtoms class.
        """
        if self.__thermo:
            self.__structures[(atoms.strainType,round(atoms.eta,3),round(atoms.V0,3))]= atoms
        elif atoms.strainType==None:
            self.__structures[(atoms.strainType,round(atoms.scale,3))]= atoms
        else:
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
    
    def set_code(self, code):
        if code in ['vasp','exciting','espresso','wien']:
            self.__code = code
        else:
            print "Unknown code '%s'. Please choose either espresso, exciting, wien or vasp"%code
            
    def get_code(self):
        return self.__code
    
    
    def status(self):
        state = Check()
        state.workdir = self.__workdir
        self.__status = state.check_calc()
        return self.__status
        
    executable    = property( fget = get_executable        , fset = set_executable)
    workdir = property( fget = get_workdir        , fset = set_workdir)
    code = property( fget = get_code       , fset = set_code)
    
    
    
    