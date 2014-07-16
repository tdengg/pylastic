Introduction
------------

**pylastic** is a python package for calculating elastic constants from first principles. Futhermore it is intended to serve as a fully automatic tool for elastic constants calculation.


Using **pylastic**' in GUI mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For automated calculations of elastic constants within a standard framework.

.. figure:: GUI_setup.png
    :width: 300px
    :align: center
    :alt: alternate text
    :figclass: align-center

Using **pylastic** as python module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Scripting with python allows for flexible usage of **pylastic**'s modules, classes and functions.

Importing vasp POSCAR files
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two options are avalilable for importing atoms positions from vasp POSCAR file:

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
	ElAtoms.poscarToAtoms(poscar)

Distorting the atoms object
^^^^^^^^^^^^^^^^^^^^^^^^^^^



The structures class
^^^^^^^^^^^^^^^^^^^^



Creating distortions
^^^^^^^^^^^^^^^^^^^^