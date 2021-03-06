"""
"""

from pylastic.analyze import Energy, Stress

#from pylastic.get_DFTdata import VASP
import pylastic.io.vasp as vasp
import pylastic.io.espresso as espresso
import pylastic.io.wien as wien
import pylastic.io.exciting as exciting
import pylastic.io.emto as emto

from pylastic.status import Check
#from pylastic.prettyPrint import FileStructure

import os
import math
import pickle
import json
import numpy as np
try:
    import matplotlib.pyplot as plt
    mpl=True
except:
    mpl=False


class ECs(Check, Energy, Stress):
    """Calculation of elastic constants and post-processing."""
    
    def __init__(self, cod='vasp', thermo=False, thmod = 'structures_phonons'):
        if cod == 'emto':
            pars = json.load(open('setup.json'))
            self.__pname = pars['emto']['pnames']['kfcd']
            self.__pname = pars['emto']['pnames']['kfcd']
            self.__emtoout = pars['emto']['jobnames']['system']+'.prn'
        #super(Check, self).__init__()
        #super(FileStructure, self).__init__()
        super(ECs, self).__init__()
        
        self.__CVS = {}
        self.__V0 = None
        self.__structures = None
        self.__cod = cod
        self.__mthd = 'Energy'
        self.__fitorder = [6,6,6]
        self.__etacalc = None
        self.__rms = {}
        self.__workdir = ''
        self.__thermodyn = thermo
        self.__mod = thmod
        self.__T = 0
        
        self.__funct = 'GGA'
        self.delPoints=False
        
        #%%%%%%%%--- CONSTANTS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        _e     = 1.602176565e-19              # elementary charge
        Bohr   = 5.291772086e-11              # a.u. to meter
        Ryd2eV = 13.605698066                 # Ryd to eV
        Angstroem = 1.e-10
        self.__cnvrtr = (_e)/(1e9*Angstroem**3.)    # Ryd/[a.u.^3] to GPa
        self.__ToGPa = (_e*Ryd2eV)/(1e9*Bohr**3.)
        self.__vToGPa = (_e)/(1e9*Angstroem**3.)
        #self.__ToGPa  = (_e*Ryd2eV)/(1e9*Bohr**3.)
        #--------------------------------------------------------------------------------------------------------------------------------
        if self.__cod=='vasp': self.__CONV=2.*self.__ToGPa
        elif self.__cod in ['espresso','exciting','wien']: self.__CONV=1.*self.__ToGPa
        
    def set_gsenergy(self, gsenergy=None):
        """Read groundstate energy for all atoms in structures.
        
        :param list gsenergy: Groundstate energy for all distorted structures (list of floats).
        
        """
        self.status()
        if not gsenergy:
            if self.__cod == 'vasp': 
                #getData = VASP()
                getData = vasp.Energy()
                outfile = 'vasprun.xml'
            elif self.__cod == 'espresso':
                getData = espresso.Energy()
                outfile = 'espresso.out'
            elif self.__cod == 'wien':
                getData = wien.Energy()
                
            elif self.__cod == 'exciting':
                getData = exciting.Energy()
                outfile = 'INFO.OUT'
            elif self.__cod == 'emto':
                getData = emto.Energy(funct=self.__funct)
                outfile = '%s/prn/%s'%(self.__pname,self.__emtoout)
            gsenergy=[]    
            for atoms in sorted(self.__structures.items()):
                
                if self.__cod == 'wien': 
                    outfile = atoms[1].path.split('/')[-1] + '.scf'
                    
                if not atoms[1].status:
                    #print atoms[1].status
                    atoms[1].gsenergy = 0
                    continue
                if atoms[1].exclude:
                    atoms[1].gsenergy_ignored = getData.get_gsenergy()
                    atoms[1].gsenergy = 0
                    continue
                if os.path.exists(atoms[1].path+'/exclude'):
                    atoms[1].gsenergy_ignored = getData.get_gsenergy()
                    atoms[1].gsenergy = 0
                    continue
                
                
                #getData.set_outfile('%s/%s/'%atoms[0] + outfile)
                #getData.set_gsEnergy()
                #print atoms[1].path, self.__workdir + '%s/%s'%(atoms[1].path.split('/')[-2],atoms[1].path.split('/')[-1])+'/' + outfile
                #getData.set_fname(self.__workdir + '%s/'%atoms[1].path.lstrip('.') + outfile)
                if 'eta' in atoms[1].path.split('/')[-1] and self.__thermodyn:getData.set_fname(self.__workdir + '%s/%s/%s'%(atoms[1].path.split('/')[-3],atoms[1].path.split('/')[-2],atoms[1].path.split('/')[-1])+'/' + outfile)
                elif 'eta' in atoms[1].path.split('/')[-1] and not self.__thermodyn:getData.set_fname(self.__workdir + '%s/%s'%(atoms[1].path.split('/')[-2],atoms[1].path.split('/')[-1])+'/' + outfile)
                else: getData.set_fname(self.__workdir + '%s'%(atoms[1].path.split('/')[-1])+'/' + outfile)
                print getData.get_fname()
                getData.set_gsenergy()
                if self.__thermodyn and self.__mod!='structures_phonons':
                    outfile_ph = 'F_TV'
                    #getData.set_fname(self.__workdir + '%s/'%atoms[1].path.lstrip('.') + outfile_ph)
                    #getData.T = self.__T
                    
                    getData.set_phenergy(self.__workdir + '%s/'%atoms[1].path.lstrip('.') + outfile_ph)
                    atoms[1].phenergy = getData.get_phenergy()
                    atoms[1].T = getData.T
                #atoms[1].gsenergy = getData.get_gsEnergy()
                    atoms[1].gsenergy = getData.get_gsenergy()/125.
                else:
                    atoms[1].gsenergy = getData.get_gsenergy()
                    gsenergy.append(atoms[1].gsenergy)
            
            if self.delPoints:
                for atoms in sorted(self.__structures.items()):
                
                    #print [atoms[1].eta for atoms in sorted(self.__structures.items())], gsenergy
                    coeff = np.polyfit([atoms[1].eta for atoms in self.__structures.items()], gsenergy, 2)
                    p = np.poly1d(coeff)
                    k=0
                    for (etas,energy) in zip(self.__structures.items(),gsenergy):
                        #print (energy-p(etas[1].eta))**2.
                        if (energy-p(etas[1].eta))**2. > 0.0004: 
                            gsenergy[k]=0.
                            atoms[1].gsenergy = 0. 
                        k+=1
        self.__gsenergy = gsenergy

    def get_gsenergy(self):
        return self.__gsenergy
    
    def set_fenergy(self, fenergy):
        
        self.__fenergy = fenergy
    
    def get_fenergy(self):
        return self.__fenergy
    
    def set_T(self,T):
        self.__T = T
    
    def get_T(self):
        return self.__T
    
    def set_stress(self, stress=None):
        """Read stress for all atoms in structures.
        :param list stress: Physical stress as list of floats
        Parameters
        ----------
            stress : 
                Physical stress.
        """
        self.status()
        if not stress:
            if self.__cod == 'vasp': 
                #getData = VASP()
                getData = vasp.Stress()
                outfile = 'vasprun.xml'
            elif self.__cod == 'espresso':
                getData = espresso.Stress()
                outfile = 'espresso.out'
            for atoms in self.__structures.items():
                
                if not atoms[1].status:
                    atoms[1].stress = np.zeros((3,3))
                    continue
                #getData.set_outfile('%s/%s/'%atoms[0] + outfile)
                #getData.set_gsEnergy()
                getData.set_fname(self.__workdir + '%s/'%atoms[1].path.lstrip('.') + outfile)
                getData.set_stress()
                #atoms[1].gsenergy = getData.get_gsEnergy()
                atoms[1].stress = getData.get_stress()
                
        self.__stress = stress
    
    def get_stress(self):
        return self.__stress
        
    def set_structures(self, excluded=None):
        """Import structures object from file."""
        self.__excluded = excluded
        if not self.__structures: 
            with open('structures.pkl', 'rb') as input:
                structures = pickle.load(input)
            if type(structures)==dict: self.__structures = structures
            else: self.__structures = structures.get_structures()
            for atoms in sorted(self.__structures.items()): atoms[1].exclude=False
        else: pass
        if self.__excluded != None: 
            for exc in self.__excluded: self.__structures[exc].exclude = True 
    def get_structures(self):
        if self.__structures == None:
            self.set_structures()
        else:
            return self.__structures
    
    def get_atomsByStraintype(self, strainType):
        """Sort structures by strain type to get an ordered list of distortions.
        
        :param str strainType: Strain type.
                
        """
        atomslist = []
        for dic in sorted(self.__structures):
            if dic[0] == strainType: 
                atomslist.append(self.__structures[dic])
        return atomslist
    
    def physicalToLagrangian(self,pstress, def_matrix):
        dm  = def_matrix
        idm = np.linalg.inv(dm)
        lstress = list(np.linalg.det(dm)*np.dot(idm,np.dot(pstress,idm)))
        return lstress
    
    def set_analytics(self):
        """Standard analysis including:
        
            * Calculation of elastic constants, 
            * Root mean square error (RMS),
            * cross validation score and
            * plot of energy vs. strain curves.
        """
        if self.__thermodyn:
            self.__A2 = {}
            
            for t in range(len(self.__T)):
                
                self.__A2[str(self.__T[t])]=[]
                self.__CVS[str(self.__T[t])]={}
                self.__rms[str(self.__T[t])]={}
                #self.status()
                # Check status of every atoms object
                strainList= self.__structures.items()[0][1].strainList
                if not len(strainList)==len(self.__fitorder): self.__fitorder=[self.__fitorder[0] for i in range(len(strainList))]
                n=1
                for stype in strainList:
                    
                    fitorder = self.__fitorder[n-1]
                    atoms = self.get_atomsByStraintype(stype)
                    self.__V0 = atoms[0].V0
                    strainList = atoms[0].strainList
                    
                    strain = [i.eta for i in atoms]
                    
                    if self.__mthd == 'Energy':
                        
                        energy = [i.gsenergy+i.phenergy[t] for i in atoms]
                        phenergy = [i.phenergy[t] for i in atoms]
                        ans = Energy(code=self.__cod)
                        ans.energy = energy
                        #if t==10 or t==100: plt.plot(strain, phenergy)
                        
                    elif self.__mthd == 'Stress':
                        
                        stress = [self.physicalToLagrangian(i.stress,i.defMatrix) for i in atoms]
                        #print stress
                        if mpl: plt.plot(strain,[i.stress[0][2] for i in atoms])
                        ans = Stress(code=self.__cod)
                        ans.set_stress(stress)
                        
                    ans.set_strain(strain)
                    
                    ans.V0 = self.__V0
                    
                    ans.set_2nd(fitorder)
                    
                    ans.set_cvs(fitorder)
                    self.__rms[str(self.__T[t])]=(ans.get_r())
                    self.__CVS[str(self.__T[t])][stype]=(ans.get_cvs())
                    self.__A2[str(self.__T[t])].append(ans.get_2nd())
                    
                    #spl = str(len(strainList))+'1'+str(n)
                    #plt.subplot(int(spl))
                    #if self.__mthd == 'energy': ans.plot_energy(fitorder=self.__fitorder)
                    n+=1
                
            if not self.__etacalc: self.__etacalc = str(strain[-1])
            if self.__mthd == 'Energy': self.set_ec(self.__etacalc)
            elif self.__mthd == 'Stress': 
                self.__sigma = ans.get_sigma()
                self.set_ec((self.__etacalc))
            
                
                
        else: 
            self.__A2 = []
            
            #self.status()
            # Check status of every atoms object
            strainList= self.__structures.items()[0][1].strainList
            n=1
            if not len(self.__fitorder)==len(strainList):
                self.__fitorder = [self.__fitorder[0] for i in range(len(strainList))]
            for stype in strainList:
                atoms = self.get_atomsByStraintype(stype)
                self.__V0 = atoms[0].V0
                strainList = atoms[0].strainList
                
                strain = [i.eta for i in atoms]
                
                if self.__mthd == 'Energy':
                    if self.__thermodyn:
                        energy = [i.gsenergy+i.phenergy for i in atoms]
                    else:    
                        energy = [i.gsenergy for i in atoms]
                    ans = Energy(code=self.__cod)
                    ans.energy = energy
                    
                elif self.__mthd == 'Stress':
                    
                    stress = [self.physicalToLagrangian(i.stress,i.defMatrix) for i in atoms]
                    
                    if mpl: plt.plot(strain,[i.stress[0][2] for i in atoms])
                    ans = Stress(code=self.__cod)
                    ans.set_stress(stress)
                    
                ans.set_strain(strain)
                
                ans.V0 = self.__V0
                
                ans.set_2nd(self.__fitorder[n-1])
                
                ans.set_cvs(self.__fitorder[n-1])
                self.__rms[stype]=ans.get_r()
                self.__CVS[stype]=ans.get_cvs()
                self.__A2.append(ans.get_2nd())
                
                #spl = str(len(strainList))+'1'+str(n)
                #plt.subplot(int(spl))
                #if self.__mthd == 'energy': ans.plot_energy(fitorder=self.__fitorder)
                n+=1
            if not self.__etacalc: self.__etacalc = str(strain[-1])
            if self.__mthd == 'Energy': self.set_ec(self.__etacalc)
            elif self.__mthd == 'Stress': 
                self.__sigma = ans.get_sigma()
                self.set_ec((self.__etacalc))
            #if mpl: plt.show()
    
    def get_rms(self):
        return self.__rms
    
    def get_CVS(self):
        return self.__CVS
    
    def plot_cvs(self, mod='F'):
        """Returns matplotlib axis instance of cross validation score plot."""
        if not mpl: raise "Problem with matplotib: Plotting not possible."
        f = plt.figure(figsize=(5,4), dpi=100)
        
        CVS = []
        strainList= self.__structures.items()[0][1].strainList
        
        if len(strainList)<=5:
            kk=1
            ll=len(strainList)
            grid=[ll]
        elif len(strainList)%5 == 0:
            kk=len(strainList)/5
            ll=5
            grid=[5 for i in range(kk)]
        else:
            kk=len(strainList)/5+1
            ll=5
            grid=[5 for i in range(kk)]
            grid[-1]=len(strainList)%5
        
        
        
        n=1
        m=1
        for stype in strainList:
            atoms = self.get_atomsByStraintype(stype)
            self.__V0 = atoms[0].V0
            strainList = atoms[0].strainList
            if self.__thermodyn and mod=='F':
                energy = [i.gsenergy+i.phenergy[-1] for i in atoms]
            elif self.__thermodyn and mod=='E0':
                energy = [i.gsenergy for i in atoms]
            elif self.__thermodyn and mod=='Fvib':
                energy = [i.phenergy[-1] for i in atoms]
            else:
                energy = [i.gsenergy for i in atoms]
            strain = [i.eta for i in atoms]
            
            spl = '1'+str(len(strainList))+str(n)
            #plt.subplot(int(spl))
            #a = f.add_subplot(int(spl))
            if (n-1)%5==0: m=0
            
            
            a = plt.subplot2grid((kk,ll), ((n-1)/5,m), colspan=1)
            
            
            j = 1
            for i in [2,4,6]:
                ans = Energy()
                ans.energy = energy
                ans.strain = strain
                ans.V0 = self.__V0
                
                fitorder = i
                ans.set_cvs(fitorder)
                CVS.append(ans.get_cvs())
                
                a.plot([cvs[1] for cvs in CVS[(n-1)*3+j-1]],[cvs[0] for cvs in CVS[(n-1)*3+j-1]], label=str(fitorder))
                a.set_title('Distortion %s'%stype)
                a.set_xlabel('strain')
                a.set_ylabel('CVS    in eV')
                a.set_xticks([0.01,0.03,0.05])
                
                j+=1
            
            n+=1
            m+=1
            
        a.legend(title='Order of fit')
        
        return f
        
    def plot_2nd(self, mod = 'F'):
        """Returns matplotlib axis instance of d2E/d(eta) plot."""
        if not mpl: raise "Problem with matplotib: Plotting not possible."
        f = plt.figure(figsize=(5,4), dpi=100)
        
        A2 = []
        
        strainList= self.__structures.items()[0][1].strainList
        
        if len(strainList)<=5:
            kk=1
            ll=len(strainList)
            grid=[ll]
        elif len(strainList)%5 == 0:
            kk=len(strainList)/5
            ll=5
            grid=[5 for i in range(kk)]
        else:
            kk=len(strainList)/5+1
            ll=5
            grid=[5 for i in range(kk)]
            grid[-1]=len(strainList)%5
        
        
        n=1
        m=1
        for stype in strainList:
            atoms = self.get_atomsByStraintype(stype)
            self.__V0 = atoms[0].V0
            strainList = atoms[0].strainList
            if self.__thermodyn and mod == 'F':
                energy = [i.gsenergy+i.phenergy[-1] for i in atoms]
            elif self.__thermodyn and mod=='E0':
                energy = [i.gsenergy for i in atoms]
            elif self.__thermodyn and mod=='Fvib':
                energy = [i.phenergy[-1] for i in atoms]
            else:
                energy = [i.gsenergy for i in atoms]
            
            strain = [i.eta for i in atoms]
            
            spl = '1'+str(len(strainList))+str(n)
            #plt.subplot(int(spl))
            #a = f.add_subplot(int(spl))
            if (n-1)%5==0: m=0
            
            
            a = plt.subplot2grid((kk,ll), ((n-1)/5,m), colspan=1)
            #print (kk,ll), ((n-1)/5,m)
            j = 0
            for i in [2,4,6]:
                ans = Energy()
                ans.energy = energy
                ans.strain = strain
                ans.V0 = self.__V0
                
                fitorder = i
                ans.set_2nd(fitorder)
                A2.append(ans.get_2nd())
                
                strains = sorted(map(float,A2[j+3*(n-1)].keys()))
                
                try:
                    dE = [A2[j+3*(n-1)][str(s)] for s in strains]
                except:
                    continue
                a.plot(strains, dE, label=str(fitorder))
                a.set_title(stype)
                a.set_xlabel('strain')
                a.set_ylabel(r'$\frac{d^2E}{d\epsilon^2}$    in eV')
                
                j+=1
            
            n+=1
            m+=1
            
        a.legend(title='Order of fit')
        return f
        
    def plot_energy(self, color=['r','g','b','c','m','y','k'], mod = 'E0'):
        """Returns matplotlib axis instance for energy-strain curve.  """
        if not mpl: raise "Problem with matplotib: Plotting not possible."
        f = plt.figure(figsize=(5,4), dpi=100)
        a = f.add_subplot(111)
        strainList= self.__structures.items()[0][1].strainList
        
        if len(strainList)<=5:
            kk=1
            ll=len(strainList)
            grid=[ll]
        elif len(strainList)%5 == 0:
            kk=len(strainList)/5
            ll=5
            grid=[5 for i in range(kk)]
        else:
            kk=len(strainList)/5+1
            ll=5
            grid=[5 for i in range(kk)]
            grid[-1]=len(strainList)%5
        
        
        n=1
        m=1
        j=0
        for stype in strainList:
            
            spl = '1'+str(len(strainList))+str(n)
            if (n-1)%5==0: m=0
            a = plt.subplot2grid((kk,ll), ((n-1)/5,m), colspan=1)
            
            
            
            fi=open(stype+'.energy','w')
            
            #self.search_for_failed()
            atoms = self.get_atomsByStraintype(stype)
            if self.__thermodyn and mod=='F':
                energy = [i.gsenergy+i.phenergy[100] for i in atoms]
            elif self.__thermodyn and mod=='E0':
                energy = [i.gsenergy for i in atoms]
            elif self.__thermodyn and mod=='Fvib':
                energy = [i.phenergy[100] for i in atoms]
            else:
                energy = [i.gsenergy for i in atoms]
            strain = [i.eta for i in atoms]
            
            ii=0
            for (e,s) in zip(energy,strain):
                if e==0.: 
                    energy.pop(ii); strain.pop(ii)
                    ii-=1
                ii+=1
            #print stype, energy, [i.scale for i in atoms]
            plt.plot(strain, energy, '%s*'%color[j%7])
            
            k=0
            for st in strain:
                fi.write('%s %s \n'%(st,energy[k]))
                k+=1
            fi.close()
            
            poly = np.poly1d(np.polyfit(strain,energy,self.__fitorder[j]))
            xp = np.linspace(min(strain), max(strain), 100)
            a.plot(xp, poly(xp),color[j%7],label=stype)
            
            a.set_title(stype)
            
            j+=1
            
            n+=1
            m+=1
        
        a.set_xlabel('strain')
        a.set_ylabel(r'energy    in eV')
        #a.legend(title='Strain type:')
        
        return f
            
    def set_fitorder(self, fitorder):
        """Set fitorder of polynomial energy - strain fit.
        
        :param list fitorder: Order of polynomial energy-strain fit (one for each distortion type).
                
        """
        self.__fitorder = fitorder
    
    def get_etacalc(self):
        return self.__etacalc
    
    def set_etacalc(self, etacalc):
        """Set maximum lagrangian strain for the fitting procedure.
        
        :param list etacalc: List of strings containing maximum lagrangian strain used for fitting.
        
        """
        self.__etacalc = etacalc
    
    def get_fitorder(self):
        return self.__fitorder
    
    def set_funct(self, funct):
        """Set exchange-correlation functional for EMTO calculations.
        
        :param str funct: Choose one of LDA, GGA or PBEsol.
                
        """
        if funct in ['LDA','GGA','PBEsol']:
            self.__funct = funct
        else:
            self.__funct = 'GGA'
            print "Unknown type of functional! Please set xc functional to one either 'LDA','GGA' or 'PBEsol'. (set to default --> GGA)"
    
    def get_funct(self):
        return self.__funct
    
    def set_ec(self, etacalc):                                           
        """Evaluate elastic constants from polynomial fitting.
        
        :param etacalc: Either list containing maximum Lagrangian strain used for fitting or only single strain used for each distortion type.
        
        """
        if not self.__thermodyn:
            C = np.zeros((6,6))
            
            LC = self.__structures.items()[0][1].LC
            if self.__mthd == 'Energy':
                if type(etacalc)==list:
                    A2=[]
                    for i in range(len(etacalc)):
                        A2.append(self.__A2[i][etacalc[i]])
                else:
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
                
            elif self.__mthd == 'Stress':
                
                if (LC == 'CI' or \
                    LC == 'CII'):
                    Matrix = np.mat([[1.0,  5.0,  0.0],
                                  [2.0,  4.0,  0.0],
                                  [3.0,  3.0,  0.0],
                                  [0.0,  0.0,  4.0],
                                  [0.0,  0.0,  5.0],
                                  [0.0,  0.0,  6.0]])
    
                if (LC == 'HI' or \
                    LC == 'HII'):
                    Matrix = np.mat([[ 1, 2, 3, 0, 0],
                                  [ 2, 1, 3, 0, 0],
                                  [ 0, 0, 3, 3, 0],
                                  [ 0, 0, 0, 0, 4],
                                  [ 0, 0, 0, 0, 5],
                                  [ 3,-3, 0, 0, 0],
                                  [ 3,-5,-1, 0, 0],
                                  [-5, 3,-1, 0, 0],
                                  [ 0, 0,-2,-1, 0],
                                  [ 0, 0, 0, 0, 6],
                                  [ 0, 0, 0, 0, 2],
                                  [-2, 2, 0, 0, 0]])
                
                if (LC == 'RI'):
                    Matrix = np.mat([[ 1, 2, 3, 4, 0, 0],
                                  [ 2, 1, 3,-4, 0, 0],
                                  [ 0, 0, 3, 0, 3, 0],
                                  [ 0, 0, 0,-1, 0, 4],
                                  [ 0, 0, 0, 6, 0, 5],
                                  [ 3,-3, 0, 5, 0, 0],
                                  [ 3,-5,-1, 6, 0, 0],
                                  [-5, 3,-1,-6, 0, 0],
                                  [ 0, 0,-2, 0,-1, 0],
                                  [ 0, 0, 0, 8, 0, 6],
                                  [ 0, 0, 0,-4, 0, 2],
                                  [-2, 2, 0, 2, 0, 0]])
                
                if (LC == 'RII'):
                    Matrix = np.mat([[ 1, 2, 3, 4, 5, 0, 0],
                                  [ 2, 1, 3,-4,-5, 0, 0],
                                  [ 0, 0, 3, 0, 0, 3, 0],
                                  [ 0, 0, 0,-1,-6, 0, 4],
                                  [ 0, 0, 0, 6,-1, 0, 5],
                                  [ 3,-3, 0, 5,-4, 0, 0],
                                  [ 3,-5,-1, 6, 2, 0, 0],
                                  [-5, 3,-1,-6,-2, 0, 0],
                                  [ 0, 0,-2, 0, 0,-1, 0],
                                  [ 0, 0, 0, 8, 4, 0, 6],
                                  [ 0, 0, 0,-4, 8, 0, 2],
                                  [-2, 2, 0, 2,-6, 0, 0]])
                
                if (LC == 'TI'):
                    Matrix = np.mat([[ 1, 2, 3, 0, 0, 0],
                                  [ 2, 1, 3, 0, 0, 0],
                                  [ 0, 0, 3, 3, 0, 0],
                                  [ 0, 0, 0, 0, 4, 0],
                                  [ 0, 0, 0, 0, 5, 0],
                                  [ 0, 0, 0, 0, 0, 6],
                                  [ 3,-5,-1, 0, 0, 0],
                                  [-5, 3,-1, 0, 0, 0],
                                  [ 0, 0,-2,-1, 0, 0],
                                  [ 0, 0, 0, 0, 6, 0],
                                  [ 0, 0, 0, 0, 2, 0],
                                  [ 0, 0, 0, 0, 0,-4]])
                
                if (LC == 'TII'):
                    Matrix = np.mat([[ 1, 2, 3, 6, 0, 0, 0],
                                  [ 2, 1, 3,-6, 0, 0, 0],
                                  [ 0, 0, 3, 0, 3, 0, 0],
                                  [ 0, 0, 0, 0, 0, 4, 0],
                                  [ 0, 0, 0, 0, 0, 5, 0],
                                  [ 0, 0, 0,-1, 0, 0, 6],
                                  [ 3,-5,-1,-4, 0, 0, 0],
                                  [-5, 3,-1, 4, 0, 0, 0],
                                  [ 0, 0,-2, 0,-1, 0, 0],
                                  [ 0, 0, 0, 0, 0, 6, 0],
                                  [ 0, 0, 0, 0, 0, 2, 0],
                                  [ 0, 0, 0, 8, 0, 0,-4]])
                
                if (LC == 'O'):
                    Matrix = np.mat([[1, 2, 3, 0, 0, 0, 0, 0, 0],
                                  [0, 1, 0, 2, 3, 0, 0, 0, 0],
                                  [0, 0, 1, 0, 2, 3, 0, 0, 0],
                                  [0, 0, 0, 0, 0, 0, 4, 0, 0],
                                  [0, 0, 0, 0, 0, 0, 0, 5, 0],
                                  [0, 0, 0, 0, 0, 0, 0, 0, 6],
                                  [3,-5,-1, 0, 0, 0, 0, 0, 0],
                                  [0, 3, 0,-5,-1, 0, 0, 0, 0],
                                  [0, 0, 3, 0,-5,-1, 0, 0, 0],
                                  [0, 0, 0, 0, 0, 0, 6, 0, 0],
                                  [0, 0, 0, 0, 0, 0, 0, 2, 0],
                                  [0, 0, 0, 0, 0, 0, 0, 0,-4],
                                  [5, 4, 6, 0, 0, 0, 0, 0, 0],
                                  [0, 5, 0, 4, 6, 0, 0, 0, 0],
                                  [0, 0, 5, 0, 4, 6, 0, 0, 0],
                                  [0, 0, 0, 0, 0, 0,-2, 0, 0],
                                  [0, 0, 0, 0, 0, 0, 0,-1, 0],
                                  [0, 0, 0, 0, 0, 0, 0, 0,-3]])
                
                if (LC == 'M'):
                    Matrix = np.mat([[ 1, 2, 3, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 1, 0, 0, 2, 3, 6, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 1, 0, 0, 2, 0, 3, 6, 0, 0, 0, 0],
                                  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 5, 0, 0],
                                  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 5, 0],
                                  [ 0, 0, 0, 1, 0, 0, 2, 0, 3, 0, 0, 0, 6],
                                  [-2, 1, 4,-5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0,-2, 0, 0, 1, 4,-5, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0,-2, 0, 0, 1, 0, 4,-5, 0, 0, 0, 0],
                                  [ 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 6, 0, 0],
                                  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 6, 0],
                                  [ 0, 0, 0,-2, 0, 0, 1, 0, 4, 0, 0,-5, 0],
                                  [ 3,-5,-1,-4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 3, 0, 0,-5,-1,-4, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 3, 0, 0,-5, 0,-1,-4, 0, 0, 0, 0],
                                  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 2, 0, 0],
                                  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 2, 0],
                                  [ 0, 0, 0, 3, 0, 0,-5, 0,-1, 0, 0,-4, 0],
                                  [-4,-6, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0,-4, 0, 0,-6, 5, 2, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0,-4, 0, 0,-6, 0, 5, 2, 0, 0, 0, 0],
                                  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-3, 0, 0],
                                  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-3, 0],
                                  [ 0, 0, 0,-4, 0, 0,-6, 0, 5, 0, 0, 2, 0],
                                  [ 5, 4, 6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 5, 0, 0, 4, 6,-3, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 5, 0, 0, 4, 0, 6,-3, 0, 0, 0, 0],
                                  [ 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0],
                                  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0],
                                  [ 0, 0, 0, 5, 0, 0, 4, 0, 6, 0, 0,-3, 0]])
                
                if (LC == 'N'):
                    Matrix = np.mat([[ 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 1, 0, 0, 0, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 4, 5, 6, 0, 0, 0],
                                  [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 4, 0, 5, 6, 0],
                                  [ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 4, 0, 5, 6],
                                  [-2, 1, 4,-3, 6,-5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0,-2, 0, 0, 0, 0, 1, 4,-3, 6,-5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 0, 4,-3, 6,-5, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0,-3, 6,-5, 0, 0, 0],
                                  [ 0, 0, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0,-3, 0, 6,-5, 0],
                                  [ 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0,-3, 0, 6,-5],
                                  [ 3,-5,-1, 6, 2,-4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 3, 0, 0, 0, 0,-5,-1, 6, 2,-4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 3, 0, 0, 0, 0,-5, 0, 0, 0,-1, 6, 2,-4, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 0, 3, 0, 0, 0, 0,-5, 0, 0, 0,-1, 0, 0, 6, 2,-4, 0, 0, 0],
                                  [ 0, 0, 0, 0, 3, 0, 0, 0, 0,-5, 0, 0, 0,-1, 0, 0, 6, 0, 2,-4, 0],
                                  [ 0, 0, 0, 0, 0, 3, 0, 0, 0, 0,-5, 0, 0, 0,-1, 0, 0, 6, 0, 2,-4],
                                  [-4,-6, 5, 1,-3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0,-4, 0, 0, 0, 0,-6, 5, 1,-3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0,-4, 0, 0, 0, 0,-6, 0, 0, 0, 5, 1,-3, 2, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 0,-4, 0, 0, 0, 0,-6, 0, 0, 0, 5, 0, 0, 1,-3, 2, 0, 0, 0],
                                  [ 0, 0, 0, 0,-4, 0, 0, 0, 0,-6, 0, 0, 0, 5, 0, 0, 1, 0,-3, 2, 0],
                                  [ 0, 0, 0, 0, 0,-4, 0, 0, 0, 0,-6, 0, 0, 0, 5, 0, 0, 1, 0,-3, 2],
                                  [ 5, 4, 6,-2,-1,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 5, 0, 0, 0, 0, 4, 6,-2,-1,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 5, 0, 0, 0, 0, 4, 0, 0, 0, 6,-2,-1,-3, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 0, 5, 0, 0, 0, 0, 4, 0, 0, 0, 6, 0, 0,-2,-1,-3, 0, 0, 0],
                                  [ 0, 0, 0, 0, 5, 0, 0, 0, 0, 4, 0, 0, 0, 6, 0, 0,-2, 0,-1,-3, 0],
                                  [ 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 4, 0, 0, 0, 6, 0, 0,-2, 0,-1,-3],
                                  [-6, 3,-2, 5,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0,-6, 0, 0, 0, 0, 3,-2, 5,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0,-6, 0, 0, 0, 0, 3, 0, 0, 0,-2, 5,-4, 1, 0, 0, 0, 0, 0, 0],
                                  [ 0, 0, 0,-6, 0, 0, 0, 0, 3, 0, 0, 0,-2, 0, 0, 5,-4, 1, 0, 0, 0],
                                  [ 0, 0, 0, 0,-6, 0, 0, 0, 0, 3, 0, 0, 0,-2, 0, 0, 5, 0,-4, 1, 0],
                                  [ 0, 0, 0, 0, 0,-6, 0, 0, 0, 0, 3, 0, 0, 0,-2, 0, 0, 5, 0,-4, 1]])
                
                sigma = np.array(self.__sigma[etacalc])
                
                ci = np.linalg.lstsq(Matrix,sigma)
                
                #-- Cubic structures ------------------------------------------------------------------------------
                if (LC == 'CI' or \
                    LC == 'CII'):
                    
                    C[0,0]=ci[0][0]
                    C[0,1]=ci[0][1]
                    C[3,3]=ci[0][2]
                    C[1,1]=C[0,0]
                    C[2,2]=C[0,0]
                    C[0,2]=C[0,1]
                    C[1,2]=C[0,1]
                    C[4,4]=C[3,3]
                    C[5,5]=C[3,3]
                
                #-- Hexagonal Structures --------------------------------------------------------------------------
                if (LC == 'HI' or \
                    LC == 'HII'):
                    C[0,0]=ci[0][0]
                    C[0,1]=ci[0][1]
                    C[0,2]=ci[0][2]
                    C[2,2]=ci[0][3]
                    C[3,3]=ci[0][4]
                    C[1,1]=C[0,0]
                    C[1,2]=C[0,2]
                    C[4,4]=C[3,3]
                    C[5,5]=0.5*(C[0,0]-C[0,1])
                
                #-- Rhombohedral I Structures ---------------------------------------------------------------------
                if (LC == 'RI'):
                    C[0,0]= ci[0][0]
                    C[0,1]= ci[0][1]
                    C[0,2]= ci[0][2]
                    C[0,3]= ci[0][3]
                    C[2,2]= ci[0][4]
                    C[3,3]= ci[0][5]
                    C[1,1]= C[0,0]
                    C[1,2]= C[0,2]
                    C[1,3]=-C[0,3]
                    C[4,5]= C[0,3]
                    C[4,4]= C[3,3]
                    C[5,5]=0.5*(C[0,0]-C[0,1])
                
                #-- Rhombohedral II Structures --------------------------------------------------------------------
                if (LC == 'RII'):
                    C[0,0]= ci[0][0]
                    C[0,1]= ci[0][1]
                    C[0,2]= ci[0][2]
                    C[0,3]= ci[0][3]
                    C[0,4]= ci[0][4]
                    C[2,2]= ci[0][5]
                    C[3,3]= ci[0][6]
                    C[1,1]= C[0,0]
                    C[1,2]= C[0,2]
                    C[1,3]=-C[0,3]
                    C[4,5]= C[0,3]
                    C[1,4]=-C[0,4]
                    C[3,5]=-C[0,4]
                    C[4,4]= C[3,3]
                    C[5,5]=0.5*(C[0,0]-C[0,1])
                
                #-- Tetragonal I Structures -----------------------------------------------------------------------
                if (LC == 'TI'):
                    C[0,0]= ci[0][0]
                    C[0,1]= ci[0][1]
                    C[0,2]= ci[0][2]
                    C[2,2]= ci[0][3]
                    C[3,3]= ci[0][4]
                    C[5,5]= ci[0][5]
                    C[1,1]= C[0,0]
                    C[1,2]= C[0,2]
                    C[4,4]= C[3,3]
                
                #-- Tetragonal II Structures ----------------------------------------------------------------------
                if (LC == 'TII'):
                    C[0,0]= ci[0][0]
                    C[0,1]= ci[0][1]
                    C[0,2]= ci[0][2]
                    C[0,5]= ci[0][3]
                    C[2,2]= ci[0][4]
                    C[3,3]= ci[0][5]
                    C[5,5]= ci[0][6]
                    C[1,1]= C[0,0]
                    C[1,2]= C[0,2]
                    C[1,5]=-C[0,5]
                    C[4,4]= C[3,3]
                
                #-- Orthorhombic Structures -----------------------------------------------------------------------
                if (LC == 'O'):
                    C[0,0]=ci[0][0]
                    C[0,1]=ci[0][1]
                    C[0,2]=ci[0][2]
                    C[1,1]=ci[0][3]
                    C[1,2]=ci[0][4]
                    C[2,2]=ci[0][5]
                    C[3,3]=ci[0][6]
                    C[4,4]=ci[0][7]
                    C[5,5]=ci[0][8]
                
                #-- Monoclinic Structures -------------------------------------------------------------------------
                if (LC == 'M'):
                    C[0,0]=ci[0][0]
                    C[0,1]=ci[0][1]
                    C[0,2]=ci[0][2]
                    C[0,5]=ci[0][3]
                    C[1,1]=ci[0][4]
                    C[1,2]=ci[0][5]
                    C[1,5]=ci[0][6]
                    C[2,2]=ci[0][7]
                    C[2,5]=ci[0][8]
                    C[3,3]=ci[0][9]
                    C[3,4]=ci[0][10]
                    C[4,4]=ci[0][11]
                    C[5,5]=ci[0][12]
                
                #-- Triclinic Structures --------------------------------------------------------------------------
                if (LC == 'N'):
                    C[0,0]=ci[0][0]
                    C[0,1]=ci[0][1]
                    C[0,2]=ci[0][2]
                    C[0,3]=ci[0][3]
                    C[0,4]=ci[0][4]
                    C[0,5]=ci[0][5]
                    C[1,1]=ci[0][6]
                    C[1,2]=ci[0][7]
                    C[1,3]=ci[0][8]
                    C[1,4]=ci[0][9]
                    C[1,5]=ci[0][10]
                    C[2,2]=ci[0][11]
                    C[2,3]=ci[0][12]
                    C[2,4]=ci[0][13]
                    C[2,5]=ci[0][14]
                    C[3,3]=ci[0][15]
                    C[3,4]=ci[0][16]
                    C[3,5]=ci[0][17]
                    C[4,4]=ci[0][18]
                    C[4,5]=ci[0][19]
                    C[5,5]=ci[0][20]
                #--------------------------------------------------------------------------------------------------
    
    
    
            for i in range(5):
                for j in range(i+1,6):
                    C[j,i] = C[i,j] 
            #%%%--- Calculating the elastic moduli ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if self.__cod in ['espresso']: C = -C/10.
            elif self.__cod  in ['vasp','exciting','wien']: C = C*self.__vToGPa/self.__V0#C = C/4.#C=C*self.__CONV/self.__V0
            elif self.__cod  in ['emto']: C = C*self.__ToGPa/self.__V0
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
        
        else:
            Cs = []
            for t in map(str,self.__T):#for t in range(len(self.__T)):
                C = np.zeros((6,6))
                
                LC = self.__structures.items()[0][1].LC
                if self.__mthd == 'Energy':
                    if type(etacalc)==list:
                        A2=[]
                        for i in range(len(etacalc)):
                            A2.append(self.__A2[t][i][etacalc[i]])
                    else:
                        if not etacalc in self.__A2[t][0].keys(): raise ValueError('Please coose one of %s'%(self.__A2[t][0].keys()))
                        A2 = [a2[etacalc] for a2 in self.__A2[t]]
                    
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
                    
                elif self.__mthd == 'Stress':
                    
                    if (LC == 'CI' or \
                        LC == 'CII'):
                        Matrix = np.mat([[1.0,  5.0,  0.0],
                                      [2.0,  4.0,  0.0],
                                      [3.0,  3.0,  0.0],
                                      [0.0,  0.0,  4.0],
                                      [0.0,  0.0,  5.0],
                                      [0.0,  0.0,  6.0]])
        
                    if (LC == 'HI' or \
                        LC == 'HII'):
                        Matrix = np.mat([[ 1, 2, 3, 0, 0],
                                      [ 2, 1, 3, 0, 0],
                                      [ 0, 0, 3, 3, 0],
                                      [ 0, 0, 0, 0, 4],
                                      [ 0, 0, 0, 0, 5],
                                      [ 3,-3, 0, 0, 0],
                                      [ 3,-5,-1, 0, 0],
                                      [-5, 3,-1, 0, 0],
                                      [ 0, 0,-2,-1, 0],
                                      [ 0, 0, 0, 0, 6],
                                      [ 0, 0, 0, 0, 2],
                                      [-2, 2, 0, 0, 0]])
                    
                    if (LC == 'RI'):
                        Matrix = np.mat([[ 1, 2, 3, 4, 0, 0],
                                      [ 2, 1, 3,-4, 0, 0],
                                      [ 0, 0, 3, 0, 3, 0],
                                      [ 0, 0, 0,-1, 0, 4],
                                      [ 0, 0, 0, 6, 0, 5],
                                      [ 3,-3, 0, 5, 0, 0],
                                      [ 3,-5,-1, 6, 0, 0],
                                      [-5, 3,-1,-6, 0, 0],
                                      [ 0, 0,-2, 0,-1, 0],
                                      [ 0, 0, 0, 8, 0, 6],
                                      [ 0, 0, 0,-4, 0, 2],
                                      [-2, 2, 0, 2, 0, 0]])
                    
                    if (LC == 'RII'):
                        Matrix = np.mat([[ 1, 2, 3, 4, 5, 0, 0],
                                      [ 2, 1, 3,-4,-5, 0, 0],
                                      [ 0, 0, 3, 0, 0, 3, 0],
                                      [ 0, 0, 0,-1,-6, 0, 4],
                                      [ 0, 0, 0, 6,-1, 0, 5],
                                      [ 3,-3, 0, 5,-4, 0, 0],
                                      [ 3,-5,-1, 6, 2, 0, 0],
                                      [-5, 3,-1,-6,-2, 0, 0],
                                      [ 0, 0,-2, 0, 0,-1, 0],
                                      [ 0, 0, 0, 8, 4, 0, 6],
                                      [ 0, 0, 0,-4, 8, 0, 2],
                                      [-2, 2, 0, 2,-6, 0, 0]])
                    
                    if (LC == 'TI'):
                        Matrix = np.mat([[ 1, 2, 3, 0, 0, 0],
                                      [ 2, 1, 3, 0, 0, 0],
                                      [ 0, 0, 3, 3, 0, 0],
                                      [ 0, 0, 0, 0, 4, 0],
                                      [ 0, 0, 0, 0, 5, 0],
                                      [ 0, 0, 0, 0, 0, 6],
                                      [ 3,-5,-1, 0, 0, 0],
                                      [-5, 3,-1, 0, 0, 0],
                                      [ 0, 0,-2,-1, 0, 0],
                                      [ 0, 0, 0, 0, 6, 0],
                                      [ 0, 0, 0, 0, 2, 0],
                                      [ 0, 0, 0, 0, 0,-4]])
                    
                    if (LC == 'TII'):
                        Matrix = np.mat([[ 1, 2, 3, 6, 0, 0, 0],
                                      [ 2, 1, 3,-6, 0, 0, 0],
                                      [ 0, 0, 3, 0, 3, 0, 0],
                                      [ 0, 0, 0, 0, 0, 4, 0],
                                      [ 0, 0, 0, 0, 0, 5, 0],
                                      [ 0, 0, 0,-1, 0, 0, 6],
                                      [ 3,-5,-1,-4, 0, 0, 0],
                                      [-5, 3,-1, 4, 0, 0, 0],
                                      [ 0, 0,-2, 0,-1, 0, 0],
                                      [ 0, 0, 0, 0, 0, 6, 0],
                                      [ 0, 0, 0, 0, 0, 2, 0],
                                      [ 0, 0, 0, 8, 0, 0,-4]])
                    
                    if (LC == 'O'):
                        Matrix = np.mat([[1, 2, 3, 0, 0, 0, 0, 0, 0],
                                      [0, 1, 0, 2, 3, 0, 0, 0, 0],
                                      [0, 0, 1, 0, 2, 3, 0, 0, 0],
                                      [0, 0, 0, 0, 0, 0, 4, 0, 0],
                                      [0, 0, 0, 0, 0, 0, 0, 5, 0],
                                      [0, 0, 0, 0, 0, 0, 0, 0, 6],
                                      [3,-5,-1, 0, 0, 0, 0, 0, 0],
                                      [0, 3, 0,-5,-1, 0, 0, 0, 0],
                                      [0, 0, 3, 0,-5,-1, 0, 0, 0],
                                      [0, 0, 0, 0, 0, 0, 6, 0, 0],
                                      [0, 0, 0, 0, 0, 0, 0, 2, 0],
                                      [0, 0, 0, 0, 0, 0, 0, 0,-4],
                                      [5, 4, 6, 0, 0, 0, 0, 0, 0],
                                      [0, 5, 0, 4, 6, 0, 0, 0, 0],
                                      [0, 0, 5, 0, 4, 6, 0, 0, 0],
                                      [0, 0, 0, 0, 0, 0,-2, 0, 0],
                                      [0, 0, 0, 0, 0, 0, 0,-1, 0],
                                      [0, 0, 0, 0, 0, 0, 0, 0,-3]])
                    
                    if (LC == 'M'):
                        Matrix = np.mat([[ 1, 2, 3, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0, 1, 0, 0, 2, 3, 6, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0, 1, 0, 0, 2, 0, 3, 6, 0, 0, 0, 0],
                                      [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 5, 0, 0],
                                      [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 5, 0],
                                      [ 0, 0, 0, 1, 0, 0, 2, 0, 3, 0, 0, 0, 6],
                                      [-2, 1, 4,-5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0,-2, 0, 0, 1, 4,-5, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0,-2, 0, 0, 1, 0, 4,-5, 0, 0, 0, 0],
                                      [ 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 6, 0, 0],
                                      [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 6, 0],
                                      [ 0, 0, 0,-2, 0, 0, 1, 0, 4, 0, 0,-5, 0],
                                      [ 3,-5,-1,-4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0, 3, 0, 0,-5,-1,-4, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0, 3, 0, 0,-5, 0,-1,-4, 0, 0, 0, 0],
                                      [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 2, 0, 0],
                                      [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 2, 0],
                                      [ 0, 0, 0, 3, 0, 0,-5, 0,-1, 0, 0,-4, 0],
                                      [-4,-6, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0,-4, 0, 0,-6, 5, 2, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0,-4, 0, 0,-6, 0, 5, 2, 0, 0, 0, 0],
                                      [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-3, 0, 0],
                                      [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-3, 0],
                                      [ 0, 0, 0,-4, 0, 0,-6, 0, 5, 0, 0, 2, 0],
                                      [ 5, 4, 6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0, 5, 0, 0, 4, 6,-3, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0, 5, 0, 0, 4, 0, 6,-3, 0, 0, 0, 0],
                                      [ 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0],
                                      [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0],
                                      [ 0, 0, 0, 5, 0, 0, 4, 0, 6, 0, 0,-3, 0]])
                    
                    if (LC == 'N'):
                        Matrix = np.mat([[ 1, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0, 1, 0, 0, 0, 0, 2, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 3, 4, 5, 6, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 4, 5, 6, 0, 0, 0],
                                      [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 4, 0, 5, 6, 0],
                                      [ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 4, 0, 5, 6],
                                      [-2, 1, 4,-3, 6,-5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0,-2, 0, 0, 0, 0, 1, 4,-3, 6,-5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 0, 4,-3, 6,-5, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0,-3, 6,-5, 0, 0, 0],
                                      [ 0, 0, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0,-3, 0, 6,-5, 0],
                                      [ 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0,-3, 0, 6,-5],
                                      [ 3,-5,-1, 6, 2,-4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0, 3, 0, 0, 0, 0,-5,-1, 6, 2,-4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0, 3, 0, 0, 0, 0,-5, 0, 0, 0,-1, 6, 2,-4, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0, 0, 3, 0, 0, 0, 0,-5, 0, 0, 0,-1, 0, 0, 6, 2,-4, 0, 0, 0],
                                      [ 0, 0, 0, 0, 3, 0, 0, 0, 0,-5, 0, 0, 0,-1, 0, 0, 6, 0, 2,-4, 0],
                                      [ 0, 0, 0, 0, 0, 3, 0, 0, 0, 0,-5, 0, 0, 0,-1, 0, 0, 6, 0, 2,-4],
                                      [-4,-6, 5, 1,-3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0,-4, 0, 0, 0, 0,-6, 5, 1,-3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0,-4, 0, 0, 0, 0,-6, 0, 0, 0, 5, 1,-3, 2, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0, 0,-4, 0, 0, 0, 0,-6, 0, 0, 0, 5, 0, 0, 1,-3, 2, 0, 0, 0],
                                      [ 0, 0, 0, 0,-4, 0, 0, 0, 0,-6, 0, 0, 0, 5, 0, 0, 1, 0,-3, 2, 0],
                                      [ 0, 0, 0, 0, 0,-4, 0, 0, 0, 0,-6, 0, 0, 0, 5, 0, 0, 1, 0,-3, 2],
                                      [ 5, 4, 6,-2,-1,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0, 5, 0, 0, 0, 0, 4, 6,-2,-1,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0, 5, 0, 0, 0, 0, 4, 0, 0, 0, 6,-2,-1,-3, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0, 0, 5, 0, 0, 0, 0, 4, 0, 0, 0, 6, 0, 0,-2,-1,-3, 0, 0, 0],
                                      [ 0, 0, 0, 0, 5, 0, 0, 0, 0, 4, 0, 0, 0, 6, 0, 0,-2, 0,-1,-3, 0],
                                      [ 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 4, 0, 0, 0, 6, 0, 0,-2, 0,-1,-3],
                                      [-6, 3,-2, 5,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0,-6, 0, 0, 0, 0, 3,-2, 5,-4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0,-6, 0, 0, 0, 0, 3, 0, 0, 0,-2, 5,-4, 1, 0, 0, 0, 0, 0, 0],
                                      [ 0, 0, 0,-6, 0, 0, 0, 0, 3, 0, 0, 0,-2, 0, 0, 5,-4, 1, 0, 0, 0],
                                      [ 0, 0, 0, 0,-6, 0, 0, 0, 0, 3, 0, 0, 0,-2, 0, 0, 5, 0,-4, 1, 0],
                                      [ 0, 0, 0, 0, 0,-6, 0, 0, 0, 0, 3, 0, 0, 0,-2, 0, 0, 5, 0,-4, 1]])
                    
                    sigma = np.array(self.__sigma[etacalc])
                    
                    ci = np.linalg.lstsq(Matrix,sigma)
                    
                    #-- Cubic structures ------------------------------------------------------------------------------
                    if (LC == 'CI' or \
                        LC == 'CII'):
                        
                        C[0,0]=ci[0][0]
                        C[0,1]=ci[0][1]
                        C[3,3]=ci[0][2]
                        C[1,1]=C[0,0]
                        C[2,2]=C[0,0]
                        C[0,2]=C[0,1]
                        C[1,2]=C[0,1]
                        C[4,4]=C[3,3]
                        C[5,5]=C[3,3]
                    
                    #-- Hexagonal Structures --------------------------------------------------------------------------
                    if (LC == 'HI' or \
                        LC == 'HII'):
                        C[0,0]=ci[0][0]
                        C[0,1]=ci[0][1]
                        C[0,2]=ci[0][2]
                        C[2,2]=ci[0][3]
                        C[3,3]=ci[0][4]
                        C[1,1]=C[0,0]
                        C[1,2]=C[0,2]
                        C[4,4]=C[3,3]
                        C[5,5]=0.5*(C[0,0]-C[0,1])
                    
                    #-- Rhombohedral I Structures ---------------------------------------------------------------------
                    if (LC == 'RI'):
                        C[0,0]= ci[0][0]
                        C[0,1]= ci[0][1]
                        C[0,2]= ci[0][2]
                        C[0,3]= ci[0][3]
                        C[2,2]= ci[0][4]
                        C[3,3]= ci[0][5]
                        C[1,1]= C[0,0]
                        C[1,2]= C[0,2]
                        C[1,3]=-C[0,3]
                        C[4,5]= C[0,3]
                        C[4,4]= C[3,3]
                        C[5,5]=0.5*(C[0,0]-C[0,1])
                    
                    #-- Rhombohedral II Structures --------------------------------------------------------------------
                    if (LC == 'RII'):
                        C[0,0]= ci[0][0]
                        C[0,1]= ci[0][1]
                        C[0,2]= ci[0][2]
                        C[0,3]= ci[0][3]
                        C[0,4]= ci[0][4]
                        C[2,2]= ci[0][5]
                        C[3,3]= ci[0][6]
                        C[1,1]= C[0,0]
                        C[1,2]= C[0,2]
                        C[1,3]=-C[0,3]
                        C[4,5]= C[0,3]
                        C[1,4]=-C[0,4]
                        C[3,5]=-C[0,4]
                        C[4,4]= C[3,3]
                        C[5,5]=0.5*(C[0,0]-C[0,1])
                    
                    #-- Tetragonal I Structures -----------------------------------------------------------------------
                    if (LC == 'TI'):
                        C[0,0]= ci[0][0]
                        C[0,1]= ci[0][1]
                        C[0,2]= ci[0][2]
                        C[2,2]= ci[0][3]
                        C[3,3]= ci[0][4]
                        C[5,5]= ci[0][5]
                        C[1,1]= C[0,0]
                        C[1,2]= C[0,2]
                        C[4,4]= C[3,3]
                    
                    #-- Tetragonal II Structures ----------------------------------------------------------------------
                    if (LC == 'TII'):
                        C[0,0]= ci[0][0]
                        C[0,1]= ci[0][1]
                        C[0,2]= ci[0][2]
                        C[0,5]= ci[0][3]
                        C[2,2]= ci[0][4]
                        C[3,3]= ci[0][5]
                        C[5,5]= ci[0][6]
                        C[1,1]= C[0,0]
                        C[1,2]= C[0,2]
                        C[1,5]=-C[0,5]
                        C[4,4]= C[3,3]
                    
                    #-- Orthorhombic Structures -----------------------------------------------------------------------
                    if (LC == 'O'):
                        C[0,0]=ci[0][0]
                        C[0,1]=ci[0][1]
                        C[0,2]=ci[0][2]
                        C[1,1]=ci[0][3]
                        C[1,2]=ci[0][4]
                        C[2,2]=ci[0][5]
                        C[3,3]=ci[0][6]
                        C[4,4]=ci[0][7]
                        C[5,5]=ci[0][8]
                    
                    #-- Monoclinic Structures -------------------------------------------------------------------------
                    if (LC == 'M'):
                        C[0,0]=ci[0][0]
                        C[0,1]=ci[0][1]
                        C[0,2]=ci[0][2]
                        C[0,5]=ci[0][3]
                        C[1,1]=ci[0][4]
                        C[1,2]=ci[0][5]
                        C[1,5]=ci[0][6]
                        C[2,2]=ci[0][7]
                        C[2,5]=ci[0][8]
                        C[3,3]=ci[0][9]
                        C[3,4]=ci[0][10]
                        C[4,4]=ci[0][11]
                        C[5,5]=ci[0][12]
                    
                    #-- Triclinic Structures --------------------------------------------------------------------------
                    if (LC == 'N'):
                        C[0,0]=ci[0][0]
                        C[0,1]=ci[0][1]
                        C[0,2]=ci[0][2]
                        C[0,3]=ci[0][3]
                        C[0,4]=ci[0][4]
                        C[0,5]=ci[0][5]
                        C[1,1]=ci[0][6]
                        C[1,2]=ci[0][7]
                        C[1,3]=ci[0][8]
                        C[1,4]=ci[0][9]
                        C[1,5]=ci[0][10]
                        C[2,2]=ci[0][11]
                        C[2,3]=ci[0][12]
                        C[2,4]=ci[0][13]
                        C[2,5]=ci[0][14]
                        C[3,3]=ci[0][15]
                        C[3,4]=ci[0][16]
                        C[3,5]=ci[0][17]
                        C[4,4]=ci[0][18]
                        C[4,5]=ci[0][19]
                        C[5,5]=ci[0][20]
                    #--------------------------------------------------------------------------------------------------
        
        
        
                for i in range(5):
                    for j in range(i+1,6):
                        C[j,i] = C[i,j] 
                #%%%--- Calculating the elastic moduli ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if self.__cod == 'espresso': C = -C/10.
                elif self.__cod in ['vasp','emto','exciting','wien']: C=C*self.__vToGPa/self.__V0#C = C/4.#C=C*self.__CONV/self.__V0#C = C/4.
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
                Cs.append(C)
            self.__C = Cs
                
    def get_ec(self):
        return self.__C
    
    def set_code(self, code):
        """Set ab-initio code.
        
        :param str code: Set ab-initio code vasp,exciting,espresso,wien or emto.
        
        """
        if code in ['vasp','exciting','espresso','wien','emto']:
            self.__cod = code
        else:
            print "Unknown code '%s'. Please choose either espresso, exciting, wien, emto or vasp"%code
            
    def get_code(self):
        return self.__cod
    
    def set_method(self, mthd):
        """Set approach of elastic constants calculation ('Energy' or 'Strain').
        
        :param str mthd: Set to 'Stress' or 'Energy'.
        
        """
        if mthd in ['Energy','Stress']: 
            self.__mthd = mthd
            
        else: print "Wrong value for method: Please choose either 'Energy' or 'Stress'!"
    
    def get_method(self):
        return self.__mthd
    
    def set_workdir(self, workdir):
        self.__workdir = workdir
    
    def get_workdir(self):
        return self.__workdir
    
    def status(self):
        state = Check()
        state.code = self.__cod
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
    ec = property( fget = get_ec        , fset = set_ec    )
    gsenergy = property( fget = get_gsenergy        , fset = set_gsenergy    )
    code = property( fget = get_code        , fset = set_code    )
    mthd = property( fget = get_method        , fset = set_method    )
    workdir = property( fget = get_workdir       , fset = set_workdir    )
    T =  property( fget = get_T       , fset = set_T    )
    functional = property( fget = get_funct       , fset = set_funct    )
            
        