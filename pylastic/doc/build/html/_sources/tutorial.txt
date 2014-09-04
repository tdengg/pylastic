Tutorials
---------

Importing vasp POSCAR files
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two options are available for importing atoms positions from vasp POSCAR file:

*	using vaspIO module (pylastic package)

	The following code sample imports the POSCAR file
	
	.. code-block:: python
	
		from pylastic.io.vasp import POS
		
		poscar = POS('POSCAR').read_pos()

		
*	using ASE (external package ASE required)

	.. code-block:: python
	
		from ase.io import vasp
		
		pos = vasp.read('POSCAR')

Creating atoms class objects from vasp POSCAR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After importing the POSCAR file one can create an ASE like atoms object which has the attributes of the POSCAR file containing information about the crystal structure (lattice vectors, number of atoms, basis, ....).

.. code-block:: python
	
	from pylastic.elatoms import ElAtoms
	from pylastic.io.vasp import POS
	
	poscar = POS('POSCAR').read_pos()
	atom = ElAtoms('vasp')
	atom.poscarToAtoms(poscar)

Distorting the atoms object
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once an atoms object is created one can make distortions by applying the *distort* method on it:

.. code-block:: python
	
	atom.distort(eta, strainType)



where the arguments `eta` and `strainType` are the lagrangian strain and the type of deformation respectively. The number of deformation types is determined by the crystal symmetry of the parent lattice.

To list all possible deformation types do

.. code-block:: python
	
	atom.strainList

The structures class
^^^^^^^^^^^^^^^^^^^^

To collect the distorted structures there is the Structures class available: 

.. code-block:: python

	from pylastic.elatoms import Structures
	
	structures = Structures('vasp')
	structures.append_structure(atom)

Starting local calculations using VASP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To start electronic structure calculations for all collected structures one has first to specify the location of the VASP binary:

.. code-block:: python

	structures.executable = '/home/t.dengg/bin/vasp/vasp.5.3/vasp'

The calculations are started after calling the method:
	
.. code-block:: python
	
	structures.calc_vasp()
	
Postprocessing
^^^^^^^^^^^^^^
The following section describes how to proceed when DFT calculations have finished. The necessary module imports for the following steps are:

.. code-block:: python

	from pylastic.elatoms import Structures, ElAtoms
	from pylastic.postprocess import ECs
	
	import matplotlib.pyplot as plt 
	
First import the structures object previously generated when setting up the calculation by calling ``set_structures``:

.. code-block:: python
	
	ec = ECs('vasp')
	ec.set_structures()
	ec.set_gsenergy()
	
With ``set_gsenergy`` the groundstate energy of each calculation is read and passed to the ``ECs`` instance as attribute.
To get 2nd order derivatives of the energy and the Cross-Validation-Score (CVS) for all distortions the ``set_analytics`` method is called: 

.. code-block:: python

	ec.set_analytics()
	
	print ec.get_CVS()
	print ec.get_rms()

For plotting matplotlib package is used. In order to plot the CVS and second energy derivative one has to call the methods ``plot_cvs`` and ``plot_2nd`` respectively. Matplotlib's ``show`` makes the figures appear on screen: 

.. code-block:: python

	import matplotlib.pyplot as plt
	
	ec.plot_cvs()
	ec.plot_2nd()
	plt.show()

Finally the elastic constants are computed for a specific maximal lagrangian strain which is given as argument and printed to standard output: 

.. code-block:: python
	
	ec.set_ec('0.05')
	print ec.get_ec()
	
	
	
	
	
	