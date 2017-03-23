Tutorials
---------

Importing vasp POSCAR files
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two options are available for importing atomic positions from *VASP* POSCAR file:

*	using ``vaspIO`` module (*pylastic* package)

	The following code imports the POSCAR file
	
	.. code-block:: python
	
		from pylastic.io.vasp import POS
		
		poscar = POS('POSCAR').read_pos()

		
*	using *ASE* (external package *ASE* required)

	.. warning:: 
	
		Not tested in new version.

	.. code-block:: python
	
		from ase.io import vasp
		
		pos = vasp.read('POSCAR')

Creating an object containing the structural information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After importing the POSCAR file one can create an object which has the attributes of the POSCAR file, containing information about the crystal structure (lattice vectors, number of atoms, basis, ....). To do so, the ``ElAtoms`` class is provided, which is similar to the *ASE* ``atoms`` class.

.. code-block:: python
	
	from pylastic.elatoms import ElAtoms
	from pylastic.io.vasp import POS
	
	poscar = POS('POSCAR').read_pos()
	atom = ElAtoms('vasp')
	atom.poscarToAtoms(poscar)

Setting calculation method
^^^^^^^^^^^^^^^^^^^^^^^^^^
	
To switch from energy (default) to the stress approach, you have to assign the ``method`` attribute to your atoms instance:

.. code-block:: python

	atom.method = 'Stress'

Distorting the atoms object
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the instance *atoms* of ``ElAtomc`` class is created, one can make distortions by applying the ``distort`` method on it:

.. code-block:: python
	
	atom.distort(eta, strainType)



where the arguments ``eta`` and ``strainType`` are the Lagrangian strain and the type of deformation respectively. The number of deformation types is determined by the crystal symmetry of the parent lattice.

To list all possible deformation types do

.. code-block:: python
	
	atom.strainList

The ``Structures`` class
^^^^^^^^^^^^^^^^^^^^^^^^

To collect the distorted structures the ``Structures`` class is available: 

.. code-block:: python

	from pylastic.elatoms import Structures
	
	structures = Structures('vasp')
	structures.append_structure(atom)

Starting local calculations using *VASP*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To start electronic structure calculations for all collected structures one has first to specify the location of the *VASP* binary:

.. code-block:: python

	structures.executable = '/home/t.dengg/bin/vasp/vasp.5.3/vasp'

The calculations are started after calling the method:
	
.. code-block:: python
	
	structures.calc_vasp()
	
Starting calculations on cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To write the paths to each singe electronic structure calculation to the file 'calcpaths':

.. code-block:: python

	f=open('calcpaths','w')
	for st in structures.get_structures().values():
		f.write(st.path.split('/')[-2]+'/'+st.path.split('/')[-1]+'/\n')
	f.close()

This file is intended for use in job subission scripts on a cluster.

	`Starting calculations on VSC3`_

Postprocessing
^^^^^^^^^^^^^^

The following section describes how to proceed when DFT calculations have finished. The modules to be imported for the subsequent steps are:

.. code-block:: python

	from pylastic.elatoms import Structures, ElAtoms
	from pylastic.postprocess import ECs
	
	import matplotlib.pyplot as plt 
	
First import the ``structures`` object previously generated when setting up the calculation by calling ``set_structures``:

.. code-block:: python
	
	ec = ECs('vasp')
	ec.set_structures()
	

.. code-block:: python	
	
	ec.set_gsenergy()
	
With ``set_gsenergy`` the groundstate energy of each calculation is read and passed to the ``ECs`` instance as attribute.
To get 2nd order derivatives of the energy and the Cross-Validation-Score (CVS) for all distortions the ``set_analytics`` method is called: 

.. code-block:: python

	ec.set_analytics()
	
	print ec.get_CVS()
	print ec.get_rms()

For plotting the package *matplotlib* is used. In order to plot the CVS and second energy derivative one has to call the methods ``plot_cvs`` and ``plot_2nd`` respectively. Matplotlib's ``show`` makes the figures appear on screen: 

.. code-block:: python

	import matplotlib.pyplot as plt
	
	ec.plot_cvs()
	ec.plot_2nd()
	plt.show()

Finally the elastic constants are computed for a specific maximal lagrangian strain which is given as argument and printed to standard output: 

.. code-block:: python
	
	ec.set_ec('0.05')
	print ec.get_ec()
	

Equation of state calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pylastic offers a tool for calculating and analyzing the equation of state (EOS) using Birch-Murnaghan as well as Vinet energy vs. volume fits. 
To generate the electronic structure input:

.. code-block:: python

	from pylastic.tools.eos import Setup, Analyze
	Setup(2.8,3.3,13, loc='cluster', cod='vasp')

After electronic structure calculations are finished you can generate an object 'analyze' containing information on the EOS

.. code-block:: python

	analyze = Analyze()

	
	
	
	