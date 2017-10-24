Introduction
------------

**pylastic** is a python package for calculating elastic constants from first principles. Furthermore it is intended to serve as a fully automatic tool for elastic constants calculation.

Concept
^^^^^^^

For calculating elastic constants two approaches are implemented:
	* Calculating 2nd derivatives of the energy with respect to strain
		.. math::

			C_{ij} = \frac{1}{V_0}\left. \frac{\partial^2 E}{\partial \eta_{i} \partial \eta_{j}} \right|_{\eta = 0} .
		
	* Derivative of stress with respect to strain
		.. math::

			C_{ij} = \frac{\partial \sigma_i}{\partial \eta_j} .
	
	For that a starting structure (parent) is distorted and energy or stress vs. strain curves are fitted with polynomials.
	
	**pylastic** allows for automation of this procedure and offers tools for data analysis such as the cross validation score (CVS).
	
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