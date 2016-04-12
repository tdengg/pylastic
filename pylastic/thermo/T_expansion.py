import numpy as np
import pickle
import json
import copy
import os

from pylastic.tools.eos import Setup, Analyze
import matplotlib.pyplot as plt

from pylastic.distort import Distort
from pylastic.elatoms import ElAtoms, Structures
from pylastic.io.vasp import POS
from pylastic.postprocess import ECs
from pylastic.tools.CVS_2nd_deriv import ANALYTICS

class Setup(Structures, Distort, POS):
    def __init__(self):
        super(Setup, self).__init__()
        self.__fname = None
        self.__poscar = None
        self.__initatom = None
    
    def deform_volume(self, Vlist, cod='vasp', executable='/home/t.dengg/bin/vasp/vasp.5.3/vasp'):
        self.__cod=cod
        if self.__cod=='vasp': from pylastic.io.vasp import POS
        


        ########################## Read in POSCAR file: ######################
        poscar = POS('POSCAR').read_pos()

        ###################### Create Structures instance: ###################
        self.__structures = Structures(self.__cod)

        ## Generate distorted structures and add them to structures object: ##
        atom = ElAtoms(self.__cod)
        atom.poscarToAtoms(poscar)
        for scale in Vlist:

            atom = ElAtoms(self.__cod)
            atom.poscarToAtoms(poscar)
            atom.scale = scale
            self.__structures.append_structure(atom)
            print atom.scale

        ####################### Write vasp input files: #######################
        self.__structures.write_structures(self.__structures)
        
    def read_POS(self, fname):
        self.__fname = fname
        self.__poscar = POS(fname).read_pos() 
        self.__initatom = ElAtoms() 
        self.__initatom.poscarToAtoms(self.__poscar)
        
    
    
    def generate_supercells(self):
        
        dirnames = [struct[1].path for struct in self.__structures.get_structures().items()]
        rootdir = os.getcwd()
        for d in dirnames:
            os.chdir(d)
            os.system('/home/MCL/t.dengg/bin/phonopy-1.7.4/bin/phonopy  -d --dim="5 5 5" -c POSCAR')
            os.system('cp POSCAR POSCAR-p')
            os.system('mv SPOSCAR POSCAR')
            os.system('pwd')
            os.chdir(rootdir)
            
