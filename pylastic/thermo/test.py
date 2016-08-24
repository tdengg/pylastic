#import imp

#vasp = imp.load_source('vasp', '/home/t.dengg/.local/lib/python2.7/site-packages/phonopy/interface/vasp.py')



import sys
print sys.path.pop(3)
sys.path.append('/home/t.dengg/.local/lib/python2.7/site-packages/phonopy')
import pickle
from phonopy import Phonopy
from phonopy.structure.atoms import Atoms as PhonopyAtoms
import numpy as np
import lxml.etree as etree
from phonopy.interface import vasp
print vasp.__file__
import os

class GET_THERMO(object):
    def __init__(self):
        
        #species = 'WRe_0.25_conv'
        species = 'WRe_0.00'
        #species = 'Re'
        ###read force constants from vasprun.xml###
        vasprun = etree.iterparse('vasprun.xml', tag='varray')
        #fc = vasp.get_force_constants_vasprun_xml(vasprun,1) #pass xml input and species atomic weight.
        ###########################################
        
        ########### read positionsl ###############
        primitive = vasp.get_atoms_from_poscar(open('POSCAR-p'),'W')
        superc =  vasp.get_atoms_from_poscar(open('POSCAR'),'W')
        ###########################################
        numbatom =  superc.get_number_of_atoms()
        #print primitive.get_cell()
        #print primitive.get_scaled_positions()
        #print superc.get_scaled_positions()
        
        print numbatom, species, os.getcwd()
        if species=='W':
        #Tungsten
            fc = vasp.get_force_constants_vasprun_xml(vasprun,1,0)
            s = 4.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['W'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
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
            
        elif species=='W_conv_2x2x2':
        #Tungsten
            fc = vasp.get_force_constants_vasprun_xml(vasprun,1,2)
            s = 2.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['W','W'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
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
        
        elif species=='W_conv':
        #Tungsten
            fc = vasp.get_force_constants_vasprun_xml(vasprun,1,2)
            s = 4.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['W','W'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
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
        
        elif species=='WRe_0.25_conv':
        #Tungsten
            fc = vasp.get_force_constants_vasprun_xml(vasprun,8,2)
            s = 4.
            a = superc.get_cell()[0][0]*2.
            print a, primitive.get_scaled_positions()
            bulk = PhonopyAtoms(symbols=['W','W'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
            phonon.set_force_constants(fc[0])
            phonon.set_dynamical_matrix()
            #print phonon.get_dynamical_matrix_at_q([0,0,0])
            mesh = [100, 100, 100]
            phonon.set_mesh(mesh)
            qpoints, weights, frequencies, eigvecs = phonon.get_mesh()
            #print frequencies
            phonon.set_total_DOS()
            
            phonon.set_thermal_properties(t_step=10,
                                          t_max=3700,
                                          t_min=0)
        
        elif species == 'WRe_B2': 
            fc = vasp.get_force_constants_vasprun_xml(vasprun,1,1)
            s = 5.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['W','Re'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
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
        elif species == 'WRe_0.00': 
            fc = vasp.get_force_constants_vasprun_xml(vasprun,1)
            s = 5.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['W'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
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
        elif species == 'WRe_0.03': 
            fc = vasp.get_force_constants_vasprun_xml(vasprun,2,0)
            s = 2.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['W'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
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
        elif species == 'WRe_0.06': 
            fc = vasp.get_force_constants_vasprun_xml(vasprun,3,0)
            s = 2.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['W'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
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
        
        elif species == 'WRe_0.09': 
            fc = vasp.get_force_constants_vasprun_xml(vasprun,4,0)
            s = 2.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['W'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
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
        elif species == 'WRe_0.12': 
            fc = vasp.get_force_constants_vasprun_xml(vasprun,5,0)
            s = 2.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['W'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
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
        
        elif species == 'WRe_0.18': 
            fc = vasp.get_force_constants_vasprun_xml(vasprun,6,0)
            s = 2.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['W'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
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
        elif species == 'WRe_0.25': 
            fc = vasp.get_force_constants_vasprun_xml(vasprun,7,0)
            s = 2.
            a = primitive.get_cell()[0][0]*2.
            print a, primitive.get_scaled_positions()
            bulk = PhonopyAtoms(symbols=['W'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
            phonon.set_force_constants(fc[0])
            phonon.set_dynamical_matrix()
            #print phonon.get_dynamical_matrix_at_q([0,0,0])
            mesh = [100, 100, 100]
            phonon.set_mesh(mesh)
            qpoints, weights, frequencies, eigvecs = phonon.get_mesh()
            print frequencies
            phonon.set_total_DOS()
            #phonon.set_partial_DOS()
            phonon.set_thermal_properties(t_step=10,
                                          t_max=3700,
                                          t_min=0)    
        
        elif species == 'WRe_0.50': 
            fc = vasp.get_force_constants_vasprun_xml(vasprun,8,0)
            s = 2.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['W'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            #print fc
            phonon.set_force_constants(fc[0])
            phonon.set_dynamical_matrix()
            #print phonon.get_dynamical_matrix_at_q([0,0,0])
            mesh = [100, 100, 100]
            phonon.set_mesh(mesh)
            qpoints, weights, frequencies, eigvecs = phonon.get_mesh()
            #print frequencies
            phonon.set_total_DOS()
            #phonon.set_partial_DOS()
            phonon.set_thermal_properties(t_step=10,
                                          t_max=3700,
                                          t_min=0)    
            
            
        elif species == 'Au': 
            fc = vasp.get_force_constants_vasprun_xml(vasprun,1,0)
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
            phonon.set_partial_DOS()
            phonon.set_thermal_properties(t_step=10,
                                          t_max=1300,
                                          t_min=0)
        
        elif species == 'Mo': 
            fc = vasp.get_force_constants_vasprun_xml(vasprun,10,0)
            s = 5.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['Mo'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
            phonon.set_force_constants(fc[0])
            phonon.set_dynamical_matrix()
            #print phonon.get_dynamical_matrix_at_q([0,0,0])
            mesh = [100, 100, 100]
            phonon.set_mesh(mesh)
            qpoints, weights, frequencies, eigvecs = phonon.get_mesh()
            print frequencies
            phonon.set_total_DOS()
            
            phonon.set_thermal_properties(t_step=10,
                                          t_max=2500,
                                          t_min=0)
        
        
        elif species == 'Re': 
            fc = vasp.get_force_constants_vasprun_xml(vasprun,9,0)
            s = 5.
            a = superc.get_cell()[0][0]*2.
            print a
            bulk = PhonopyAtoms(symbols=['Re'] * 1,
                                scaled_positions= primitive.get_scaled_positions())
            bulk.set_cell(np.diag((a, a, a)))
            phonon = Phonopy(bulk,
                             [[s,0.,0.],[0.,s,0.],[0.,0.,s]],
                             primitive_matrix=[[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],
                             distance=0.01, factor=15.633302)
            print fc
            phonon.set_force_constants(fc[0])
            phonon.set_dynamical_matrix()
            #print phonon.get_dynamical_matrix_at_q([0,0,0])
            mesh = [100, 100, 100]
            phonon.set_mesh(mesh)
            qpoints, weights, frequencies, eigvecs = phonon.get_mesh()
            print frequencies
            phonon.set_total_DOS()
            
            phonon.set_thermal_properties(t_step=10,
                                          t_max=2500,
                                          t_min=0)
        
        
        f = open('F_TV','w')
        for t, free_energy, entropy, cv in np.array(phonon.get_thermal_properties()).T:
            #print t, cv
            #print ("%12.3f " + "%15.7f" * 3) % ( t, free_energy, entropy, cv )
            f.write(("%12.3f " + "%15.7f" + "\n") % ( t, free_energy))
        f.close()
        
        fc = open('thermal_properties','w')
        for t, free_energy, entropy, cv in np.array(phonon.get_thermal_properties()).T:
            fc.write(("%12.3f " + "%15.7f" *3 + "\n") % ( t, free_energy, entropy, cv ))
        fc.close()
        
        #phonon.plot_thermal_properties().show()
        
        #phonon.plot_total_DOS().show()
        phonon.write_total_DOS()
        #phonon.write_partial_DOS()
        phonon.write_yaml_thermal_properties()
        
        bands = []
        
        #### PRIMITIVE
        
        q_start  = np.array([0.0, 0.0, 0.0])
        #q_start  = np.array([0.5, 0.5, 0.0])
        q_end    = np.array([-0.5, 0.5, 0.5])
        #q_end    = np.array([0., 0., 0.])
        band = []
        for i in range(101):
            band.append(q_start + (q_end - q_start) / 100 * i)
        bands.append(band)
        
        band = []
        
        
        q_start  = np.array([-0.5, 0.5, 0.5])
        #q_start  = np.array([0., 0., 0.])
        q_end    = np.array([0.25, 0.25, 0.25])
        #q_end    = np.array([1., 0., 0.])
        
        for i in range(101):
            #band.append([-0.5+3*1/400*i, 0.5-1/400*i, 0.5-1/400*i])
            band.append(q_start + (q_end - q_start) / 100 * i)
        bands.append(band)
        #print band
        q_start  = np.array([0.25, 0.25, 0.25])
        q_end    = np.array([0., 0., 0.])
        band = []
        for i in range(101):
            band.append(q_start + (q_end - q_start) / 100 * i)
        bands.append(band)
        
        q_start  = np.array([0., 0., 0.])
        q_end    = np.array([0.0, 0., 0.5])
        band = []
        for i in range(101):
            band.append(q_start + (q_end - q_start) / 100 * i)
        bands.append(band)
        
        #q_start  = np.array([0.0, 0.0, 0.0])
        #q_end    = np.array([0.5, 0.5, 0.5])
        #band = []
        #for i in range(101):
        #    band.append(q_start + (q_end - q_start) / 100 * i)
        #bands.append(band)
        
        """
        ###### CONVENTIONAL CELL ######
        q_start  = np.array([0.0, 0.0, 0.0])
        q_end    = np.array([-0.5, 0.5, 0.5])
        band = []
        for i in range(101):
            band.append(q_start + (q_end - q_start) / 100 * i)
        bands.append(band)
        
        
        q_start  = np.array([-0.5, 0.5, 0.5])
        q_end    = np.array([1./4., 1./4., 1./4.])
        band = []
        for i in range(101):
            band.append(q_start + (q_end - q_start) / 100 * i)
        bands.append(band)
        q_start  = np.array([1./4., 1./4., 1./4.])
        q_end    = np.array([0.0, 0.0, 0.0])
        band = []
        for i in range(101):
            band.append(q_start + (q_end - q_start) / 100 * i)
        bands.append(band)
        q_start  = np.array([0.0, 0.0, 0.0])
        q_end    = np.array([0.5, 0.0, 0.0])
        band = []
        for i in range(101):
            band.append(q_start + (q_end - q_start) / 100 * i)
        bands.append(band)
        """
        phonon.set_band_structure(bands)
        #phonon.plot_band_structure().show()
        
        q_points, distances, frequencies, eigvecs = phonon.get_band_structure()
        disp = {'q':q_points, 'distances':distances, 'frequencies':frequencies, 'eigvecs':eigvecs}
        f = open('ph_dispersion.pkl','w')
        pickle.dump(disp, f)
        f.close()