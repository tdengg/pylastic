Starting calculations on VSC3
-----------------------------

* Start with copying your calculation directory to the cluster.

* Change to the root directory of your calculations.

* Generate job submission script named 'run_vasp'

::

	#!/bin/tcsh
	#SBATCH --qos=devel_0128
	#SBATCH --partition=mem_0128
	#SBATCH -N 1 
	#SBATCH -J $SLURM_SUBMIT_DIR:t 
	#SBATCH --tasks-per-node=16
	#SBATCH -t 00:10:00
	mpirun /home/lv70372/tde/bin/vasp/vasp.5.3.3-collin_DNGZhalf_mkl_vtst_3.0d/vasp
	

* Write shell script named 'run.sh' to submit all jobs:

::

	#!/bin/bash
	#
	
	
	for path in 3*; do
		cd $path
		cdir="$(pwd)"
		
	    cd ${line}
	    cp ../run_vasp .
	    sbatch run_vasp
		
		cd ..
	done



* Start job script

::

	./run.sh