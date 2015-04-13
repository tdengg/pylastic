import sys
print sys.path.pop(2)
import pickle
from phonopy import Phonopy
from phonopy.structure.atoms import Atoms as PhonopyAtoms
import numpy as np
import lxml.etree as etree
from phonopy.interface import vasp
print vasp.__file__

class GET_THERMO(object):
    def __init__(self):
        
        ###read force constants from vasprun.xml###
        vasprun = etree.iterparse('vasprun.xml', tag='varray')
        fc = vasp.get_force_constants_vasprun_xml(vasprun,0) #pass xml input and species atomic weight.
        ###########################################
        
        ########### read positionsl ###############
        primitive = vasp.get_atoms_from_poscar(open('POSCAR-p'),'W')
        superc =  vasp.get_atoms_from_poscar(open('POSCAR'),'W')
        ###########################################
        numbatom =  superc.get_number_of_atoms()
        #print primitive.get_cell()
        #print primitive.get_scaled_positions()
        #print superc.get_scaled_positions()
        """ #Tungsten
        s = 5.
        a = superc.get_cell()[0][0]*2.
        bulk = PhonopyAtoms(symbols=['W'] * 1,
                            scaled_positions= primitive.get_scaled_positions())
        bulk.set_cell(np.diag((a, a, a)))
        phonon = Phonopy(bulk,
                         [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                         primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                         distance=0.01, factor=15.633302)
        
        phonon.set_force_constants(fc[0])
        phonon.set_dynamical_matrix()
        #print phonon.get_dynamical_matrix_at_q([0,0,0])
        mesh = [100, 100, 100]
        phonon.set_mesh(mesh)
        qpoints, weights, frequencies, eigvecs = phonon.get_mesh()
        print frequencies
        phonon.set_total_DOS()
        
        phonon.set_thermal_properties(t_step=10,
                                      t_max=3700,
                                      t_min=0)
        """
        
        #Gold
        s = 5.
        a = superc.get_cell()[0][0]
        bulk = PhonopyAtoms(symbols=['Au'] * 1,
                            scaled_positions= primitive.get_scaled_positions())
        bulk.set_cell(np.diag((a, a, a)))
        phonon = Phonopy(bulk,
                         [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                         primitive_matrix=[[0.5, 0.5, 0.0],[0.0, 0.5, 0.5],[0.5, 0.0, 0.5]],
                         distance=0.01, factor=15.633302)
        
        phonon.set_force_constants(fc[0])
        phonon.set_dynamical_matrix()
        #print phonon.get_dynamical_matrix_at_q([0,0,0])
        mesh = [100, 100, 100]
        phonon.set_mesh(mesh)
        qpoints, weights, frequencies, eigvecs = phonon.get_mesh()
        print frequencies
        phonon.set_total_DOS()
        
        phonon.set_thermal_properties(t_step=10,
                                      t_max=1300,
                                      t_min=0)
        
        f = open('F_TV','w')
        for t, free_energy, entropy, cv in np.array(phonon.get_thermal_properties()).T:
            print t, cv
            print ("%12.3f " + "%15.7f" * 3) % ( t, free_energy, entropy, cv )
            f.write(("%12.3f " + "%15.7f" + "\n") % ( t, free_energy))
        f.close()
        
        fc = open('thermal_properties','w')
        for t, free_energy, entropy, cv in np.array(phonon.get_thermal_properties()).T:
            fc.write(("%12.3f " + "%15.7f" *3 + "\n") % ( t, free_energy, entropy, cv ))
        fc.close()
        
        #phonon.plot_thermal_properties().show()
        
        #phonon.plot_total_DOS().show()
        phonon.write_total_DOS()
        
        phonon.write_yaml_thermal_properties()
        
        bands = []
        q_start  = np.array([0.0, 0.0, 0.0])
        q_end    = np.array([0.0, 0.5, 0.5])
        band = []
        for i in range(101):
            band.append(q_start + (q_end - q_start) / 100 * i)
        bands.append(band)
        
        
        q_start  = np.array([0.0, 0.5, 0.5])
        q_end    = np.array([3./8., 3./4., 3./8.])
        band = []
        for i in range(101):
            band.append(q_start + (q_end - q_start) / 100 * i)
        bands.append(band)
        q_start  = np.array([3./8., 3./4., 3./8.])
        q_end    = np.array([0.0, 0.0, 0.0])
        band = []
        for i in range(101):
            band.append(q_start + (q_end - q_start) / 100 * i)
        bands.append(band)
        q_start  = np.array([0.0, 0.0, 0.0])
        q_end    = np.array([0.5, 0.5, 0.5])
        band = []
        for i in range(101):
            band.append(q_start + (q_end - q_start) / 100 * i)
        bands.append(band)
        
        phonon.set_band_structure(bands)
        #phonon.plot_band_structure().show()
        
        q_points, distances, frequencies, eigvecs = phonon.get_band_structure()
        disp = {'q':q_points, 'distances':distances, 'frequencies':frequencies, 'eigvecs':eigvecs}
        f = open('ph_dispersion.pkl','w')
        pickle.dump(disp, f)
        f.close()