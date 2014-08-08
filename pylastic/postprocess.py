"""
"""

from pylastic.analyze import Energy
from pylastic.get_DFTdata import VASP
from pylastic.status import Check
from pylastic.prettyPrint import FileStructure


import json
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

class ECs_old(object):
    def __init__(self):
        #%%%%%%%%--- CONSTANTS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        _e     = 1.602176565e-19              # elementary charge
        #Bohr   = 5.291772086e-11              # a.u. to meter
        #Ryd2eV = 13.605698066                 # Ryd to eV
        Angstroem = 1.e-10
        self.__cnvrtr = (_e)/(1e9*Angstroem**3)    # Ryd/[a.u.^3] to GPa
        #--------------------------------------------------------------------------------------------------------------------------------
        self.__CVS
        self.V0 = None
        self.__cod = 'vasp'
        
    def set_CVS(self):
        """Calculate cross validation score."""
        getData = VASP()
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
            
            ans = Energy(strain,energy,self.V0)
            ans.set_2nd(6)
            print ans.get_2nd()
        
    def get_CVS(self):
        return self.__CVS
            
class ECs(Check, FileStructure):
    """Calculate elastic constants, CVS calculation, post-processing."""
    
    def __init__(self):
        
        super(Check, self).__init__()
        super(FileStructure, self).__init__()
        
        self.__CVS = []
        self.__V0 = None
        self.__structures = None
        self.__cod = 'vasp'
        self.__fitorder = 4
        self.__etacalc = '0.04'
        self.__workdir = './'
        #%%%%%%%%--- CONSTANTS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        _e     = 1.602176565e-19              # elementary charge
        #Bohr   = 5.291772086e-11              # a.u. to meter
        #Ryd2eV = 13.605698066                 # Ryd to eV
        Angstroem = 1.e-10
        self.__cnvrtr = (_e)/(1e9*Angstroem**3)    # Ryd/[a.u.^3] to GPa
        #--------------------------------------------------------------------------------------------------------------------------------
        
    def set_gsenergy(self, gsenergy=None):
        """Read groundstate energy for all atoms in structures.
        
        Parameters
        ----------
            gsenergy : float???? DEBUG???? 
                Groundstate energy.
        """
        self.status()
        if not gsenergy:
            getData = VASP()
            for atoms in self.__structures.items():
                # Do interpolation if calculation failed...... COMPLETE IT
                
                if not atoms[1].status:
                    atoms[1].gsenergy = 0
                    continue
                getData.set_outfile('%s/%s/vasprun.xml'%atoms[0])
                getData.set_gsEnergy()
                atoms[1].gsenergy = getData.get_gsEnergy()
                
        
                
                
        self.__gsenergy = gsenergy
    
    def get_gsenergy(self):
        return self.__gsenergy
    
        
    def set_structures(self):
        """Import structures object from file."""
        if not self.__structures: 
            with open('structures.pkl', 'rb') as input:
                self.__structures = pickle.load(input).get_structures()
        else: pass
        
    def get_structures(self):
        if not self.__structures:
            self.set_structures()
        else:
            return self.__structures
    
    def get_atomsByStraintype(self, strainType):
        """Sort structures by strain type to get an ordered list of distortions.
        
        Parameters
        ----------
            strainType : string 
                Strain type.
        """
        atomslist = []
        for dic in sorted(self.__structures):
            if dic[0] == strainType: 
                atomslist.append(self.__structures[dic])
        return atomslist
    
    def set_analytics(self):
        """Standard analysis including:
            * Calculation of elastic constants, 
            * cross validation score and 
            * plot of the energy strain curves.
        """
        self.__A2 = []
        
        #self.status()
        # Check status of every atoms object
        strainList= self.__structures.items()[0][1].strainList
        n=1
        for stype in strainList:
            atoms = self.get_atomsByStraintype(stype)
            self.__V0 = atoms[0].V0
            strainList = atoms[0].strainList
            energy = [i.gsenergy for i in atoms]
            strain = [i.eta for i in atoms]
            
            ans = Energy(strain,energy,self.__V0)
            ans.set_2nd(self.__fitorder)
            ans.set_cvs(self.__fitorder)
            self.__CVS.append(ans.get_cvs())
            self.__A2.append(ans.get_2nd())
            
            #spl = str(len(strainList))+'1'+str(n)
            #plt.subplot(int(spl))
            #ans.plot_energy()
            n+=1
        
        self.set_C(self.__etacalc)
        #plt.show()
    
    def get_CVS(self):
        return self.__CVS
        
    def plot_cvs(self):
        """Returns matplotlib axis instance of cross validation score plot."""
        f = plt.figure(figsize=(5,4), dpi=100)
        
        CVS = []
        strainList= self.__structures.items()[0][1].strainList
        n=1
        for stype in strainList:
            atoms = self.get_atomsByStraintype(stype)
            self.__V0 = atoms[0].V0
            strainList = atoms[0].strainList
            energy = [i.gsenergy for i in atoms]
            strain = [i.eta for i in atoms]
            
            spl = '1'+str(len(strainList))+str(n)
            #plt.subplot(int(spl))
            a = f.add_subplot(int(spl))
            j = 1
            for i in [2,4,6]:
                ans = Energy(strain,energy,self.__V0)
                self.__fitorder = i
                ans.set_cvs(self.__fitorder)
                CVS.append(ans.get_cvs())
                
                a.plot([cvs[1] for cvs in CVS[(n-1)*3+j-1]],[cvs[0] for cvs in CVS[(n-1)*3+j-1]], label=str(self.__fitorder))
                a.set_title(stype)
                a.set_xlabel('strain')
                a.set_ylabel('CVS')
                
                j+=1
            
            n+=1
            
        a.legend(title='Order of fit')
        
        return f
        
    def plot_2nd(self):
        """Returns matplotlib axis instance of d2E/d(eta) plot."""
        f = plt.figure(figsize=(5,4), dpi=100)
        
        A2 = []
        
        strainList= self.__structures.items()[0][1].strainList
        n=1
        for stype in strainList:
            atoms = self.get_atomsByStraintype(stype)
            self.__V0 = atoms[0].V0
            strainList = atoms[0].strainList
            energy = [i.gsenergy for i in atoms]
            strain = [i.eta for i in atoms]
            
            spl = '1'+str(len(strainList))+str(n)
            a = f.add_subplot(int(spl))
            
            j = 0
            for i in [2,4,6]:
                ans = Energy(strain,energy,self.__V0)
                self.__fitorder = i
                ans.set_2nd(self.__fitorder)
                A2.append(ans.get_2nd())
                strains = sorted(map(float,A2[j].keys()))
                dE = [A2[j][str(s)] for s in strains]
                
                a.plot(strains, dE, label=str(self.__fitorder))
                a.set_title(stype)
                a.set_xlabel('strain')
                a.set_ylabel('dE')
                
                j+=1
            
            n+=1
            
        a.legend(title='Order of fit')
        return f
        
        
    def set_fitorder(self, fitorder):
        """Set fitorder of polynomial energy - strain fit.
        
        Parameters
        ----------
            fitorder : integer 
                Order of polynomial energy-strain fit.
        """
        self.__fitorder = fitorder
    
    def get_etacalc(self):
        return self.__etacalc
    
    def set_etacalc(self, etacalc):
        """Set maximum lagrangian strain for the fitting procedure.
        
        Parameters
        ----------
            etacalc : float 
                Maximum lagrangian strain.
        """
        self.__etacalc = etacalc
    
    def get_fitorder(self):
        return self.__fitorder
    
    def set_C(self, etacalc):
        """Evaluate elastic constants from polynomials.
        
        Parameters
        ----------
            etacalc : float 
                Maximum lagrangian strain.
            A2 : list
                Coefficients of polynomial fit.
        """
        C = np.zeros((6,6))
        
        LC = self.__structures.items()[0][1].LC
        if not etacalc in self.__A2[0].keys(): raise ValueError('Please coose one of %s'%(self.__A2[0].keys()))
        A2 = [a2[etacalc] for a2 in self.__A2]
        
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
        
        
        for i in range(5):
            for j in range(i+1,6):
                C[j,i] = C[i,j] 
        #%%%--- Calculating the elastic moduli ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        self.BV = (C[0,0]+C[1,1]+C[2,2]+2*(C[0,1]+C[0,2]+C[1,2]))/9
        self.GV = ((C[0,0]+C[1,1]+C[2,2])-(C[0,1]+C[0,2]+C[1,2])+3*(C[3,3]+C[4,4]+C[5,5]))/15
        self.EV = (9*self.BV*self.GV)/(3*self.BV+self.GV)
        self.nuV= (1.5*self.BV-self.GV)/(3*self.BV+self.GV)
        self.S  = np.linalg.inv(C)
        self.BR = 1/(self.S[0,0]+self.S[1,1]+self.S[2,2]+2*(self.S[0,1]+self.S[0,2]+self.S[1,2]))
        self.GR =15/(4*(self.S[0,0]+self.S[1,1]+self.S[2,2])-4*(self.S[0,1]+self.S[0,2]+self.S[1,2])+3*(self.S[3,3]+self.S[4,4]+self.S[5,5]))
        self.ER = (9*self.BR*self.GR)/(3*self.BR+self.GR)
        self.nuR= (1.5*self.BR-self.GR)/(3*self.BR+self.GR)
        self.BH = 0.50*(self.BV+self.BR)
        self.GH = 0.50*(self.GV+self.GR)
        self.EH = (9.*self.BH*self.GH)/(3.*self.BH+self.GH)
        self.nuH= (1.5*self.BH-self.GH)/(3.*self.BH+self.GH)
        self.AVR= 100.*(self.GV-self.GR)/(self.GV+self.GR)
        #--------------------------------------------------------------------------------------------------------------------------------
        self.__C = C
        
    def get_C(self):
        return self.__C
    
    def status(self):
        state = Check()
        state.workdir = self.__workdir
        state.structures = self.__structures
        self.__status, self.__structuresinst, statusstring = state.check_calc()
        try:
            self.__structures = self.__structuresinst.get_structures()
        except:
            self.__structures = self.__structuresinst
    
    structures    = property( fget = get_structures         , fset = set_structures    )
    fitorder    = property( fget = get_fitorder        , fset = set_fitorder    )
    etacalc    = property( fget = get_etacalc        , fset = set_etacalc    )
    C = property( fget = get_C        , fset = set_C    )
    gsenergy = property( fget = get_gsenergy        , fset = set_gsenergy    )
    
    
if __name__ == '__main__':
    ECs_old().set_CVS()
            
        