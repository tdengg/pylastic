Tutorials
---------

Importing vasp POSCAR files
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two options are available for importing atoms positions from vasp POSCAR file:

*	using vaspIO module (pylastic package)

	The following code sample imports the POSCAR file
	
	.. code-block:: python
	
		from pylastic.vaspIO import POS
		
		poscar = POS.read_pos('POSCAR')
		
*	using ASE (external package ASE required)

	.. code-block:: python
	
		from ase.io import vasp
		
		pos = vasp.read('POSCAR')

Creating atoms class objects from vasp POSCAR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After importing the POSCAR file one can create an ASE like atoms object which has the attributes of the POSCAR file containing information about the crystal structure (lattice vectors, number of atoms, basis, ....).

.. code-block:: python
	
	from pylastic.elatoms import ElAtoms
	from pylastic.vaspIO import POS
	
	poscar = POS.read_pos('POSCAR')
	atom = ElAtoms()
	atom.poscarToAtoms(poscar)

Distorting the atoms object
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once an atoms object is created one can make distortions by using the distort() method on it:

.. code-block:: python
	
	atom.distort(eta, strainType)

where the arguments eta and strainType are the strain and the type of deformation. The number of deformation types is determined by the crystal symmetry of the parent lattice.

To list all possible deformation types do

.. code-block:: python
	
	atom.strainList

The structures class
^^^^^^^^^^^^^^^^^^^^

To collect the distorted structures there is the Structures class available: 

.. code-block:: python

	from pylastic.elatoms import Structures
	
	structures = Structures()
	structures.append_structure(atom)
