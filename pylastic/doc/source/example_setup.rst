Example: 2nd order elastic constants at the groundstate.
--------

Setup:
______

.. code-block:: python

	import numpy as np
	from elatoms import Structures, ElAtoms
	from vaspIO import POS
	
	########################## Read in POSCAR file: ######################
	poscar = POS('POSCAR').read_pos()
	
	###################### Create Structures instance: ###################
	structures = Structures()
	
	## Generate distorted structures and add them to structures object: ##
	atom = ElAtoms()
	atom.poscarToAtoms(poscar)
	for etan in np.linspace(-0.05,0.05,11):
		
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

Postprocessing:
_______________

.. code-block:: python

	from elatoms import Structures, ElAtoms
	from postprocess import ECs
	
	ec = ECs()
	ec.set_structures()
	ec.set_gsenergy()
	ec.set_analytics()