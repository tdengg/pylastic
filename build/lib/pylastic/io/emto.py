import os
import re
import sys
import numpy as np

class POS(object):
    def __init__(self, verbouse=False):
        self.__verbouse=verbouse
        self.__basevect=[]
        self.__lattvect_pos=None
        self.__natoms = 0
        self.__pos={}
        self.__NQ3=1
        self.__N_species+=0
    def read_in(self):
        
        return
    
    def read_kstr(self, fname):
        if fname: self.__kstr = open(fname)
        
        lines = self.__kstr.readlines()
        self.__split_kstr =  [re.split(r'(\s+)', l) for l in lines]
        i=0
        j=0
        k=0
        for line in self.__split_kstr:
            if line[0].startswith('BS'): 
                self.__pos['vlatt_%s'%i] = [line[2],line[6],line[10]]
                if self.__lattvect_pos==None: self.__lattvect_pos=j
                i+=1
            elif line[0].startswith('QX'):
                self.__natoms+=1
                self.__basevect.append([line[2],line[6],line[10]])
                k+=1
            elif line[0].startswith('NQ3'):
                self.__NQ3=int(line[2])
            j+=1
        
        return
    
    def write_kstr(self):
        
        return
    
    def read_shape(self, fname):
        if fname: self.__shape = open(fname)
        lines = self.__shape.readlines()
        self.__split_shape =  [re.split(r'(\s+)', l) for l in lines]
        return
    
    def write_shape(self):
        
        return
    
    def read_kgrn(self, fname):
        if fname: self.__kgrn = open(fname)
        lines = self.__kgrn.readlines()
        self.__split_kgrn =  [re.split(r'(\s+)', l) for l in lines]
        i=0
        j=0
        k=0
        for line in self.__split_kstr:
            if line[0].startswith('SWS'): 
                self.__SWS = line[2]
            elif line[0].startswith('Symb'):
                i+=1
            elif 0<i<(self.__NQ3+1): 
                self.__N_species+=1
                i+=1
            
            
            
            j+=1
        return
    
    def write_kgrn(self):
        
        return
    
    def read_kfcd(self, fname):
        
        if fname: self.__kfcd = open(fname)
        lines = self.__kfcd.readlines()
        self.__split_kfcd =  [re.split(r'(\s+)', l) for l in lines]
        return
    
    def write_sgroup(self, pos):
        
        f = open('sgroup.in','w')
        
        f.write('P\n')
        
        a = np.sqrt((pos["vlatt_1"][0])**2. + (pos["vlatt_1"][1])**2. + (pos["vlatt_1"][2])**2.) * pos["scale"]
        b = np.sqrt((pos["vlatt_2"][0])**2. + (pos["vlatt_2"][1])**2. + (pos["vlatt_2"][2])**2.) * pos["scale"]
        c = np.sqrt((pos["vlatt_3"][0])**2. + (pos["vlatt_3"][1])**2. + (pos["vlatt_3"][2])**2.) * pos["scale"]
        alpha = np.arccos(np.dot(pos["vlatt_2"]*pos["scale"],pos["vlatt_3"]*pos["scale"])/(b*c))*180./np.pi
        beta = np.arccos(np.dot(pos["vlatt_3"]*pos["scale"],pos["vlatt_1"]*pos["scale"])/(a*c))*180./np.pi
        gamma = np.arccos(np.dot(pos["vlatt_1"]*pos["scale"],pos["vlatt_2"]*pos["scale"])/(a*b))*180./np.pi
        
        f.write(str(a) +' '+ str(b) +' '+ str(c) +' '+ str(alpha) +' '+ str(beta) +' '+ str(gamma) + '\n')
        natom = 0
        for s in pos['natoms']: natom += s
        f.write(str(natom)+'\n')
        if pos["csystem"] in ['d','D']:
            for i in range(len(pos["natoms"])):
                for j in range(pos["natoms"][i]):
                    
                    f.write(pos["vbasis"]["species_" + str(i+1)][j][0]  + pos["vbasis"]["species_" + str(i+1)][j][1]  + pos["vbasis"]["species_" + str(i+1)][j][2] + '\n')
                    f.write('Species_' + str(i+1) + '\n')
        else: print 'Basis vectors in Cartesian coordinates not supported yet!!! \n NOT WRITTEN TO sgroup.in!!!!!'
        f.close()
        
class Energy(object):
    """Get energies from EMTO output file."""
    def __init__(self, fname = 'vasprun.xml'):
        self.__fname = fname
        
    def set_T(self, T):
        self.__T = T
        
    def get_T(self):
        return self.__T

    def set_gsenergy(self):
        """Get groundstate energy from vasprun.xml"""
        self.__gsenergy = float(0)
        
        
    def get_gsenergy(self):
        """Return groundstate energy."""
        return self.__gsenergy
    
    def set_phenergy(self, ph_file):
        """Read phonon free energy from phonopy output (default filename: F_TV) for temperature T."""
        g = open(ph_file)
        lines = g.readlines()
        phenergy = []
        T = []
        for temp in lines:
            phenergy.append(float(temp.split()[1])/96.47244)
            T.append(float(temp.split()[0]))
        g.close()
        
        self.__phenergy = phenergy
        self.__T=T
    
    def get_phenergy(self):
        """Return phonon free energy for temperature T"""
        return self.__phenergy
    
    def set_fname(self,fname):
        self.__fname = fname
        
    def get_fname(self):
        return self.__fname
    
    fname = property( fget = get_fname        , fset = set_fname)
    T = property( fget = get_T       , fset = set_T)
    
    