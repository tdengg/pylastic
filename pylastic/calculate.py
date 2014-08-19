import numpy as np
from elatoms import Structures, ElAtoms
from pylastic.io.vasp import POS

NoP = 11
etamax = 0.05

########################## Read in POSCAR file: ######################
poscar = POS('POSCAR').read_pos()

###################### Create Structures instance: ###################
structures = Structures()

## Generate distorted structures and add them to structures object: ##
atom = ElAtoms()
atom.poscarToAtoms(poscar)
for etan in np.linspace(-etamax,etamax,NoP):
    
    for strains in range(len(atom.strainList)):
        atom = ElAtoms()
        atom.poscarToAtoms(poscar)
        atom.distort(eta=etan, strainType_index = strains)
        structures.append_structure(atom)
        
####################### Write vasp input files: #######################
structures.write_structures()

#################### Start local vasp calculation: ####################
structures.executable = '/home/t.dengg/bin/vasp/vasp.5.3/vasp'
structures.calc_vasp()
