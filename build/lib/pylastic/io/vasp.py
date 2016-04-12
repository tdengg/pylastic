import numpy as np
import lxml.etree as et
import os

class POS(object):
    def __init__(self, fname='POSCAR'):
        self.__fname = fname
        if fname: self.car = open(fname)
        self.__verbose = True
        self.__T = 0
        
    def read_pos(self):
        if self.__verbose: print "Reading input file: '%s' ...."%(self.__fname)
        p_dict = {}
        p_dict["name"] = self.car.readline()
        if self.__fname: p_dict["path"] = self.__fname
        else: p_dict["path"] = os.getcwd()+'/POSCAR'
        
        scale = self.lta()[0]
        
        try:
            scale = np.float(scale)
        except:
            print 'Found list of lattice parameters'
        p_dict["scale"] = scale #Volume if < 0
        
        p_dict["vlatt_1"] = np.array(map(np.float,self.lta()))
        p_dict["vlatt_2"] = np.array(map(np.float,self.lta()))
        p_dict["vlatt_3"] = np.array(map(np.float,self.lta()))
        p_dict["natoms"] = map(int,self.lta())
        
        selective = self.car.readline()
        if selective.split()[0][0] in ['s','S']: 
            p_dict["selective"] = True
            p_dict["csystem"] = self.car.readline().split()[0][0]
        else: 
            p_dict["selective"] = False
            p_dict["csystem"] = selective.split()[0][0]
        p_dict["vbasis"] = {}
        n=1
        for i in range(len(p_dict["natoms"])):
            #p_dict["vbasis"]["b_"+str(i+1)] = []
            for j in range(p_dict["natoms"][i]):
                p_dict["vbasis"]["b_"+str(n)]=np.array(map(np.float,self.lta()))
                n+=1
        if self.__verbose: print "'%s' read in as dictionary."%(self.__fname) 
        return p_dict
        #self.write_sgroup(p_dict)
    
    def read_in(self):
        i_dict = {}
        newlines = []
        lines = self.car.readlines()
        for l in lines: newlines.append(l.strip().split(' = '))
        #print newlines
        for pars in newlines:
            if not pars == ['']: i_dict[pars[0]] = pars[1]
        
        #print i_dict
        return i_dict
        
        
    def lta(self):
        line = self.car.readline()
        l = line.split()
        
        for i in range(len(l)):
            if l[i] in [' ','\n','']: l.pop(i)
            else: continue
        
        return l
    
    def write_pos(self, pos, fileName):
        if self.__verbose: print "Writing file '%s'"%fileName
        posout = open(fileName, 'w')
        posout.write('COMMENT!\n')
        posout.write(str(pos['scale']) + '\n')
        for i in range(3): posout.write(str(pos['vlatt_1'][i]) + ' ')
        posout.write('\n')
        for i in range(3): posout.write(str(pos['vlatt_2'][i]) + ' ')
        posout.write('\n')
        for i in range(3): posout.write(str(pos['vlatt_3'][i]) + ' ')
        posout.write('\n')
        for n in pos['natoms']: posout.write(str(n) + ' ')
        posout.write('\n')
        
        if pos['selective']: posout.write('Selective\n')
        posout.write(pos['csystem'] + '\n')
        
        for n in range(len(pos["vbasis"])): 
            
            for s in pos['vbasis']["b_%s"%(str(n+1))]:
                
                #for l in s: 
                posout.write(str(s)+' ')#!!!!!!!
            posout.write('\n')
                
        posout.close()
        
    def write_in(self, inpars, fileName):
        inout = open(fileName, 'w')
        for par in inpars: 
            inout.write(par + ' = ' + str(inpars[par]) + '\n')
        
        inout.close()
        
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
        
        for i in range(len(pos["natoms"])):
            for j in range(pos["natoms"][i]):
                if i>0:
                    b0_i=pos["vbasis"]["b_" + str(i)]
                    b1_i=pos["vbasis"]["b_" + str(i+1)]
                    if b0_i[0]==b1_i[0]and b0_i[1]==b1_i[1]and b0_i[2]==b1_i[2]:
                        natom = natom-1
                        continue
        
        f.write(str(natom)+'\n')
        
        if pos["csystem"] in ['d','D']:
            n=1
            for i in range(len(pos["natoms"])):
                for j in range(pos["natoms"][i]):
                    if i>0:
                        b0_i=pos["vbasis"]["b_" + str(n-1)]
                        b1_i=pos["vbasis"]["b_" + str(n)]
                        if b0_i[0]==b1_i[0]and b0_i[1]==b1_i[1]and b0_i[2]==b1_i[2]:
                            natom = natom-1
                            continue
                    
                    f.write(str(pos["vbasis"]["b_" + str(n)][0])  + ' '  + str(pos["vbasis"]["b_" + str(n)][1])  + ' ' + str(pos["vbasis"]["b_" + str(n)][2]) + '\n')
                    f.write('Species_' + str(i+1) + '\n')
                    n+=1
        else: print 'Basis vectors in Cartesian coordinates not supported yet!!! \n NOT WRITTEN TO sgroup.in!!!!!'
        f.close()
    
    def set_fname(self, fname):
        self.__fname = fname
    
    def get_fname(self):
        return self.__fname
    
    fname = property( fget = get_fname        , fset = set_fname)

        
    
#    def write_sgroup_ase(self, atoms):
#        
#        f = open('sgroup.in','w')
#        
#        f.write('P\n')
#        
#        a = np.sqrt((atoms.cell[0][0])**2. + (atoms.cell[0][1])**2. + (atoms.cell[0][2])**2.) * pos["scale"]
#        b = np.sqrt((atoms.cell[1][0])**2. + (atoms.cell[1][1])**2. + (atoms.cell[1][2])**2.) * pos["scale"]
#        c = np.sqrt((atoms.cell[2][0])**2. + (atoms.cell[2][1])**2. + (atoms.cell[2][2])**2.) * pos["scale"]
#        alpha = np.arccos(np.dot(atoms.cell[1]*pos["scale"],atoms.cell[2]*pos["scale"])/(b*c))*180./np.pi
#        beta = np.arccos(np.dot(atoms.cell[2]*pos["scale"],atoms.cell[0]*pos["scale"])/(a*c))*180./np.pi
#        gamma = np.arccos(np.dot(atoms.cell[0]*pos["scale"],atoms.cell[1]*pos["scale"])/(a*b))*180./np.pi
#        
#        f.write(str(a) +' '+ str(b) +' '+ str(c) +' '+ str(alpha) +' '+ str(beta) +' '+ str(gamma) + '\n')
#        natom = 0
#        for s in pos['natoms']: natom += s
#        f.write(str(natom)+'\n')
#        if pos["csystem"] in ['d','D']:
#            for i in range(len(pos["natoms"])):
#                for j in range(pos["natoms"][i]):
#                    
#                    f.write(pos["vbasis"]["species_" + str(i+1)][j][0] + ' ' + pos["vbasis"]["species_" + str(i+1)][j][1] + ' ' + pos["vbasis"]["species_" + str(i+1)][j][2] + '\n')
#                    f.write('Species_' + str(i+1) + '\n')
#        else: print 'Basis vectors in Cartesian coordinates not supported yet!!! \n NOT WRITTEN TO sgroup.in!!!!!'
#        f.close()


class Energy():
    """Get energies from espresso output file."""
    def __init__(self, fname = 'vasprun.xml'):
        self.__fname = fname
        
    def set_T(self, T):
        self.__T = T
        
    def get_T(self):
        return self.__T

    def set_gsenergy(self):
        """Get groundstate energy from vasprun.xml"""
        vasprun = et.parse(self.__fname)
        elem = vasprun.xpath("//scstep/energy[i/@name='hartreedc']/i[@name = 'e_0_energy']")
        self.__gsenergy = float(elem[1].text)
        
        
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
    
class Stress():
    def __init__(self, fname = 'vasprun.xml'):
        self.__fname = fname
        self.__stress = None
        
    def set_stress(self):
        if (os.path.exists(self.__fname)):
            stress_m = []
          
            vaspout = et.parse(self.__fname)
            stress_in = vaspout.xpath("//varray[@name = 'stress']/v")
            for stress_i in stress_in: stress_m.append(map(float, (stress_i.text).split()))
            self.__stress = stress_m
                   
        
    def get_stress(self):
        return self.__stress
    
    def set_fname(self,fname):
        self.__fname = fname
        
    def get_fname(self):
        return self.__fname
    
    fname = property( fget = get_fname        , fset = set_fname)
    
    
if __name__ == "__main__": POS('INCAR').read_in()