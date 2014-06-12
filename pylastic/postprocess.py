import analyze
import get_DFTdata
import json
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

class ECs_old(object):
    def __init__(self):
        #%%%%%%%%--- CONSTANTS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        _e     = 1.602176565e-19              # elementary charge
        #Bohr   = 5.291772086e-11              # a.u. to meter
        #Ryd2eV = 13.605698066                 # Ryd to eV
        Angstroem = 1.e-10
        self.__cnvrtr = (_e)/(1e9*Angstroem**3)    # Ryd/[a.u.^3] to GPa
        #--------------------------------------------------------------------------------------------------------------------------------
        
        self.V0 = None
        self.__cod = 'vasp'
        
    def set_CVS(self):
        
        getData = get_DFTdata.VASP()
        f=open('info.json')
        dic = json.load(f)
        
        for key in sorted(dic.keys()):
            energy = []
            strain = []
            for key2 in sorted(map(int,dic[key].keys())):
                
                getData.set_outfile('Dst%.2d'%int(key)+'/Dst%.2d'%int(key)+'_%.2d'%int(key2)+'/vasprun.xml')
                getData.set_gsEnergy()
                energy.append(getData.get_gsEnergy())
                strain.append(dic[key][str(key2)]['eta'])
            self.V0 = dic[key][str(key2)]['V0']
            
            ans = analyze.Energy(strain,energy,self.V0)
            ans.set_2nd(6)
            print ans.get_2nd()
            
class ECs(object):
    def __init__(self):
        self.__V0 = None
        self.__structures = None
        self.__cod = 'vasp'
        self.__fitorder = 6
        self.__etacalc = '0.03'
        #%%%%%%%%--- CONSTANTS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        _e     = 1.602176565e-19              # elementary charge
        #Bohr   = 5.291772086e-11              # a.u. to meter
        #Ryd2eV = 13.605698066                 # Ryd to eV
        Angstroem = 1.e-10
        self.__cnvrtr = (_e)/(1e9*Angstroem**3)    # Ryd/[a.u.^3] to GPa
        #--------------------------------------------------------------------------------------------------------------------------------
        
    def set_gsenergy(self, gsenergy=None):
        if not gsenergy:
            getData = get_DFTdata.VASP()
            for atoms in self.__structures.items():
                getData.set_outfile('%s/%s/vasprun.xml'%atoms[0])
                getData.set_gsEnergy()
                atoms[1].gsenergy = getData.get_gsEnergy()
                
                
        self.__gsenergy = gsenergy
    
    def get_gsenergy(self):
        return self.__gsenergy
    
        
    def set_structures(self):
        if not self.__structures: 
            with open('structures.pkl', 'rb') as input:
                self.__structures = pickle.load(input)
        else: pass
        
    def get_structures(self):
        if not self.__structures:
            self.set_structures()
        else:
            return self.__structures
    
    def get_atomsByStraintype(self, strainType):
        atomslist = []
        for dic in sorted(self.__structures):
            if dic[0] == strainType: 
                atomslist.append(self.__structures[dic])
        return atomslist
    
    def set_analytics(self):
        A2 = []
        
        strainList= self.__structures.items()[0][1].strainList
        n=1
        for stype in strainList:
            atoms = self.get_atomsByStraintype(stype)
            self.__V0 = atoms[0].V0
            strainList = atoms[0].strainList
            energy = [i.gsenergy for i in atoms]
            strain = [i.eta for i in atoms]
            
            ans = analyze.Energy(strain,energy,self.__V0)
            ans.set_2nd(self.__fitorder)
            ans.set_cvs(self.__fitorder)
            A2.append(ans.get_2nd())
            
            spl = str(len(strainList))+'1'+str(n)
            plt.subplot(int(spl))
            ans.plot_energy()
            n+=1
        print A2
        self.set_ECs(A2, self.__etacalc)
        plt.show()
        
    def set_fitorder(self, fitorder):
        self.__fitorder = fitorder
    
    def get_etacalc(self):
        return self.__etacalc
    
    def set_etacalc(self, etacalc):
        self.__etacalc = etacalc
    
    def get_fitorder(self):
        return self.__fitorder
    
    def set_ECs(self, A2, etacalc):
        C = np.zeros((6,6))
        
        LC = self.__structures.items()[0][1].LC
        
        A2 = [a2[etacalc] for a2 in A2]
        print A2
        #%%%--- Cubic structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (LC == 'CI' or \
            LC == 'CII'):
            C[0,0] =-2.*(A2[0]-3.*A2[1])/3.
            C[1,1] = C[0,0]
            C[2,2] = C[0,0]
            C[3,3] = A2[2]/6.
            C[4,4] = C[3,3]
            C[5,5] = C[3,3]
            C[0,1] = (2.*A2[0]-3.*A2[1])/3.
            C[0,2] = C[0,1]
            C[1,2] = C[0,1]
        #--------------------------------------------------------------------------------------------------------------------------------
        
        #%%%--- Hexagonal structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (LC == 'HI' or \
            LC == 'HII'):
            C[0,0] = 2.*A2[3]
            C[0,1] = 2./3.*A2[0] + 4./3.*A2[1] - 2.*A2[2] - 2.*A2[3]
            C[0,2] = 1./6.*A2[0] - 2./3.*A2[1] + 0.5*A2[2]
            C[1,1] = C[0,0]
            C[1,2] = C[0,2]
            C[2,2] = 2.*A2[2]
            C[3,3] =-0.5*A2[2] + 0.5*A2[4]
            C[4,4] = C[3,3]
            C[5,5] = .5*(C[0,0] - C[0,1])
        #--------------------------------------------------------------------------------------------------------------------------------
        
        #%%%--- Rhombohedral I structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (LC == 'RI'):
            C[0,0] = 2.*A2[3]
            C[0,1] = A2[1]- 2.*A2[3]
            C[0,2] = .5*( A2[0] - A2[1] - A2[2])
            C[0,3] = .5*(-A2[3] - A2[4] + A2[5])
            C[1,1] = C[0,0]
            C[1,2] = C[0,2]
            C[1,3] =-C[0,3]
            C[2,2] = 2.*A2[2]
            C[3,3] = .5*A2[4]
            C[4,4] = C[3,3]
            C[4,5] = C[0,3]
            C[5,5] = .5*(C[0,0] - C[0,1])
        #--------------------------------------------------------------------------------------------------------------------------------
        
        #%%%--- Rhombohedral II structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (LC == 'RII'):
            C[0,0] = 2.*A2[3]
            C[0,1] = A2[1]- 2.*A2[3]
            C[0,2] = .5*( A2[0] - A2[1] - A2[2])
            C[0,3] = .5*(-A2[3] - A2[4] + A2[5])
            C[0,4] = .5*(-A2[3] - A2[4] + A2[6])
            C[1,1] = C[0,0]
            C[1,2] = C[0,2]
            C[1,3] =-C[0,3]
            C[1,4] =-C[0,4]    
            C[2,2] = 2.*A2[2]
            C[3,3] = .5*A2[4]
            C[3,5] =-C[0,4]
            C[4,4] = C[3,3]
            C[4,5] = C[0,3]
            C[5,5] = .5*(C[0,0] - C[0,1])
        #--------------------------------------------------------------------------------------------------------------------------------
        
        #%%%--- Tetragonal I structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (LC == 'TI'):
            C[0,0] = (A2[0]+2.*A2[1])/3.+.5*A2[2]-A2[3]
            C[0,1] = (A2[0]+2.*A2[1])/3.-.5*A2[2]-A2[3]
            C[0,2] = A2[0]/6.-2.*A2[1]/3.+.5*A2[3]
            C[1,1] = C[0,0]
            C[1,2] = C[0,2]
            C[2,2] = 2.*A2[3]
            C[3,3] = .5*A2[4]
            C[4,4] = C[3,3]
            C[5,5] = .5*A2[5]
        #--------------------------------------------------------------------------------------------------------------------------------
        
        #%%%--- Tetragonal II structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (LC == 'TII'):
            C[0,0] = (A2[0]+2.*A2[1])/3.+.5*A2[2]-A2[4]
            C[1,1] = C[0,0]
            C[0,1] = (A2[0]+2.*A2[1])/3.-.5*A2[2]-A2[4]
            C[0,2] = A2[0]/6.-(2./3.)*A2[1]+.5*A2[4]
            C[0,5] = (-A2[2]+A2[3]-A2[6])/4.
            C[1,2] = C[0,2]
            C[1,5] =-C[0,5]
            C[2,2] = 2.*A2[4]
            C[3,3] = .5*A2[5]
            C[4,4] = C[3,3]
            C[5,5] = .5*A2[6]
        #--------------------------------------------------------------------------------------------------------------------------------
        
        #%%%--- Orthorhombic structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (LC == 'O'):
            C[0,0] = 2.*A2[0]/3.+4.*A2[1]/3.+A2[3]-2.*A2[4]-2.*A2[5]
            C[0,1] = 1.*A2[0]/3.+2.*A2[1]/3.-.5*A2[3]-A2[5]
            C[0,2] = 1.*A2[0]/3.-2.*A2[1]/3.+4.*A2[2]/3.-.5*A2[3]-A2[4]
            C[1,1] = 2.*A2[4]
            C[1,2] =-2.*A2[1]/3.-4.*A2[2]/3.+.5*A2[3]+A2[4]+A2[5]
            C[2,2] = 2.*A2[5]
            C[3,3] = .5*A2[6]
            C[4,4] = .5*A2[7]
            C[5,5] = .5*A2[8]
        #--------------------------------------------------------------------------------------------------------------------------------
        
        #%%%--- Monoclinic structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (LC == 'M'):
            C[0,0] = 2.*A2[0]/3.+8.*(A2[1]+A2[2])/3.-2.*(A2[5]+A2[8]+A2[9])
            C[0,1] = A2[0]/3.+4.*(A2[1]+A2[2])/3.-2.*A2[5]-A2[9]
            C[0,2] =(A2[0]-4.*A2[2])/3.+A2[5]-A2[8]
            C[0,5] =-1.*A2[0]/6.-2.*(A2[1]+A2[2])/3.+.5*(A2[5]+A2[7]+A2[8]+A2[9]-A2[12])
            C[1,1] = 2.*A2[8]
            C[1,2] =-4.*(2.*A2[1]+A2[2])/3.+2.*A2[5]+A2[8]+A2[9]+A2[12]
            C[1,5] =-1.*A2[0]/6.-2.*(A2[1]+A2[2])/3.-.5*A2[3]+A2[5]+.5*(A2[7]+A2[8]+A2[9])
            C[2,2] = 2.*A2[9]
            C[2,5] =-1.*A2[0]/6.+2.*A2[1]/3.-.5*(A2[3]+A2[4]-A2[7]-A2[8]-A2[9]-A2[12])
            C[3,3] = .5*A2[10]
            C[3,4] = .25*(A2[6]-A2[10]-A2[11])
            C[4,4] = .5*A2[11]
            C[5,5] = .5*A2[12]
        #--------------------------------------------------------------------------------------------------------------------------------
        
        #%%%--- Triclinic structures ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (LC == 'N'):
            C[0,0] = 2.*A2[0]
            C[0,1] = 1.*(-A2[0]-A2[1]+A2[6])
            C[0,2] = 1.*(-A2[0]-A2[2]+A2[7])
            C[0,3] = .5*(-A2[0]-A2[3]+A2[8]) 
            C[0,4] = .5*(-A2[0]+A2[9]-A2[4])
            C[0,5] = .5*(-A2[0]+A2[10]-A2[5])
            C[1,1] = 2.*A2[1]
            C[1,2] = 1.*(A2[11]-A2[1]-A2[2])
            C[1,3] = .5*(A2[12]-A2[1]-A2[3])
            C[1,4] = .5*(A2[13]-A2[1]-A2[4])
            C[1,5] = .5*(A2[14]-A2[1]-A2[5])
            C[2,2] = 2.*A2[2] 
            C[2,3] = .5*(A2[15]-A2[2]-A2[3])
            C[2,4] = .5*(A2[16]-A2[2]-A2[4])
            C[2,5] = .5*(A2[17]-A2[2]-A2[5])
            C[3,3] = .5*A2[3]
            C[3,4] = .25*(A2[18]-A2[3]-A2[4])
            C[3,5] = .25*(A2[19]-A2[3]-A2[5])
            C[4,4] = .5*A2[4]
            C[4,5] = .25*(A2[20]-A2[4]-A2[5])
            C[5,5] = .5*A2[5]
        #--------------------------------------------------------------------------------------------------------------------------------
        
        #%%%--- Calculating the elastic moduli ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (self.__cod == 'WIEN2k'):
            CONV = self.__cnvrtr * 1.
        if (self.__cod == 'exciting'):
            CONV = self.__cnvrtr * 1.
        if (self.__cod == 'ESPRESSO'):
            CONV = self.__cnvrtr * 1.
        if (self.__cod == 'vasp') or (self.__cod == 'vasp_T'):
            CONV = self.__cnvrtr * 1.
        for i in range(5):
            for j in range(i+1,6):
                C[j,i] = C[i,j] 
        
        C = C * CONV/self.__V0
        
        BV = (C[0,0]+C[1,1]+C[2,2]+2*(C[0,1]+C[0,2]+C[1,2]))/9
        GV = ((C[0,0]+C[1,1]+C[2,2])-(C[0,1]+C[0,2]+C[1,2])+3*(C[3,3]+C[4,4]+C[5,5]))/15
        EV = (9*BV*GV)/(3*BV+GV)
        nuV= (1.5*BV-GV)/(3*BV+GV)
        S  = np.linalg.inv(C)
        BR = 1/(S[0,0]+S[1,1]+S[2,2]+2*(S[0,1]+S[0,2]+S[1,2]))
        GR =15/(4*(S[0,0]+S[1,1]+S[2,2])-4*(S[0,1]+S[0,2]+S[1,2])+3*(S[3,3]+S[4,4]+S[5,5]))
        ER = (9*BR*GR)/(3*BR+GR)
        nuR= (1.5*BR-GR)/(3*BR+GR)
        BH = 0.50*(BV+BR)
        GH = 0.50*(GV+GR)
        EH = (9.*BH*GH)/(3.*BH+GH)
        nuH= (1.5*BH-GH)/(3.*BH+GH)
        AVR= 100.*(GV-GR)/(GV+GR)
        #--------------------------------------------------------------------------------------------------------------------------------
        print C[0,0], C[0,1], C[3,3]
    
    structures    = property( fget = get_structures         , fset = set_structures    )
    
    
if __name__ == '__main__':
    ECs_old().set_CVS()
            
        