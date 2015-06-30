Example script using pylastic
--------

The following example will show the procedure of calculating elastic constants (at T=0K) using either the GUI version or a *python* script.

a.) Using the GUI
^^^^^^^^^^^^^^^^^

b.) Using python modules in a script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Setup:
______

.. code-block:: python

	import numpy as np
	from pylastic.elatoms import Structures, ElAtoms
	from pylastic.io.vasp import POS
	
	########################## Read in POSCAR file: ######################
	poscar = POS('POSCAR').read_pos()
	
	###################### Create Structures instance: ###################
	structures = Structures('vasp')
	
	## Generate distorted structures and add them to structures object: ##
	atom = ElAtoms('vasp')
	atom.poscarToAtoms(poscar)
	for etan in np.linspace(-0.05,0.05,11):
		
		for strains in range(len(atom.strainList)):
			atom = ElAtoms('vasp')
			atom.poscarToAtoms(poscar)
			atom.distort(eta=etan, strainType_index = strains)
			structures.append_structure(atom)
			
	####################### Write vasp input files: #######################
	structures.write_structures(structures)
	
	#################### Start local vasp calculation: ####################
	structures.executable = '/home/t.dengg/bin/vasp/vasp.5.3/vasp'
	structures.calc_vasp()

Postprocessing:
_______________

.. code-block:: python

	from pylastic.elatoms import Structures, ElAtoms
	from pylastic.postprocess import ECs
	
	ec = ECs('vasp')
	ec.set_structures()
	ec.set_gsenergy()
	ec.set_analytics()