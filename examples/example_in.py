import sys
sys.path.append('/home/t.dengg/git/pylastic/pylastic')

import numpy as np
from pylastic.elatoms import Structures, ElAtoms
from pylastic.io.exciting import INPUT


########################## Read in POSCAR file: ######################
poscar = INPUT('input.xml').read_in()

###################### Create Structures instance: ###################
structures = Structures()
structures.code = 'exciting'

## Generate distorted structures and add them to structures object: ##
atom = ElAtoms()
atom.code = 'exciting'
atom.poscarToAtoms(poscar)
for etan in np.linspace(-0.05,0.05,11):
	
	for strains in range(len(atom.strainList)):
		atom = ElAtoms()
		atom.code = 'exciting'
		atom.poscarToAtoms(poscar)
		atom.distort(eta=etan, strainType_index = strains)
		structures.append_structure(atom)
		
####################### Write vasp input files: #######################
structures.write_structures(structures)

#################### Start local vasp calculation: ####################
#structures.executable = '/home/t.dengg/bin/vasp/vasp.5.3/vasp'
#structures.calc_vasp()

