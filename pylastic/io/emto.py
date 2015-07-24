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
        
        self.__format = {"kgrn":{"SWS":{"read":(8,14) , "write":"{0:7.4f}"}}, "kstr":{"A.....":{"read":(10,19) , "write":"{0:10.7f}"},"B.....":{"read":(30,39) , "write":"{0:10.7f}"},"C.....":{"read":(50,59) , "write":"{0:10.7f}"},"BSX":{"read":(10,19) , "write":"{0:10.7f}"},"BSY":{"read":(30,39) , "write":"{0:10.7f}"},"BSZ":{"read":(50,59) , "write":"{0:10.7f}"}}}
        
    def find(self, ifile, fname, parname):
        j=0
        val=[]
        for line in ifile:
            if parname in line:
                f_read=self.__format[fname][parname]["read"]
                val.append(float(line[f_read[0]:f_read[1]]))
            
                j+=1
        if len(val)==1: val=val[0]
        
        return val
    
    def replace(self, ifile, fname, parname, val):
        i=0
        j=0
        for line in ifile:
            if parname in line:
        
                f_read=self.__format[fname][parname]["read"]
                li = list(line)
                del li[f_read[0]:f_read[1]]
                if type(val) is list: li[f_read[0]] = self.__format[fname][parname]["write"].format(float(val[j]))
                else: li[f_read[0]] = self.__format[fname][parname]["write"].format(float(val))
        
                line = ''.join(li)
                ifile[i] = line
            
                j+=1
            i+=1
        return ifile
    
    def read_in(self):
        
        return
    
    def read_kstr(self, fname):
        if fname: kstr = open(fname)
        
        self.__kstr = kstr.readlines()
        A = self.find(self.__kstr, "kstr", "A.....")
        B = self.find(self.__kstr, "kstr", "B.....")
        C = self.find(self.__kstr, "kstr", "C.....")
        
        BSX = self.find(self.__kstr, "kstr", "BSX")
        BSY = self.find(self.__kstr, "kstr", "BSY")
        BSZ = self.find(self.__kstr, "kstr", "BSZ")
        
        for i in range(len(BSX)):
            self.__pos['vlatt_%s'%i][0] = BSX[i]
            self.__pos['vlatt_%s'%i][1] = BSY[i]
            self.__pos['vlatt_%s'%i][2] = BSZ[i]
            
            
        return
    
    def write_kstr(self, fname):
        BSX = np.zeros(3)
        BSY = np.zeros(3)
        BSZ = np.zeros(3)
        
        for i in range(len(BSX)):
            BSX[i] = self.__pos['vlatt_%s'%i][0]
            BSY[i] = self.__pos['vlatt_%s'%i][1]
            BSZ[i] = self.__pos['vlatt_%s'%i][2]
        
        self.__kstr = self.replace(self.__kstr,"kstr","BSX", BSX)
        self.__kstr = self.replace(self.__kstr,"kstr","BSY", BSY)
        self.__kstr = self.replace(self.__kstr,"kstr","BSZ", BSZ)
        return
    
    def read_shape(self, fname):
        if fname: shape = open(fname)
        self.__shape = shape.readlines()
        
        
        return
    
    def write_shape(self):
        
        return
    
    def read_kgrn(self, fname):
        if fname: kgrn = open(fname)
        self.__kgrn = kgrn.readlines()
        
        SWS = self.find(self.__kgrn,"kgrn","SWS")
        self.__pos['scale'] = SWS
        return
    
    def write_kgrn(self):
        SWS = self.__pos['scale']
        self.__kgrn = self.replace(self.__kgrn,"kgrn","SWS", SWS)
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
    def __init__(self, fname = 'etot.dat', funct='PBEsol'):
        self.__funct = funct
        self.__fname = fname
        
    def set_T(self, T):
        self.__T = T
        
    def get_T(self):
        return self.__T

    def set_gsenergy(self):
        """Get groundstate energy from vasprun.xml"""
        if self.__funct == 'LDA': ind = 1
        elif self.__funct == 'GGA': ind = 2
        elif self.__funct == 'PBEsol': ind = 3
        for line in open(self.__fname,'r'): self.__gsenergy = float(line.split()[ind])
        
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
    
    