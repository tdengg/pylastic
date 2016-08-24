import numpy as np
import pickle
import json
import copy
import os

from pylastic.distort import Distort
from pylastic.elatoms import ElAtoms, Structures
from pylastic.io.vasp import POS
from pylastic.postprocess import ECs
from pylastic.tools.CVS_2nd_deriv import ANALYTICS

class Setup(Structures, Distort, POS):
    def __init__(self):
        super(Setup, self).__init__()
        self.__fname = None
        self.__poscar = None
        
        self.__initatom = None
        
        
    def read_POS(self, fname):
        self.__fname = fname
        self.__poscar = POS(fname).read_pos() 
        self.__initatom = ElAtoms() 
        self.__initatom.poscarToAtoms(self.__poscar)
        
    def setup_DFT(self, etamax, N):
        self.__structures = Structures(thermo=True)
        atom = copy.deepcopy(self.__initatom)
        for eta in np.linspace(-etamax,etamax,N):
            atom = copy.deepcopy(self.__initatom)
            atom.deform_volume(eta)
            #self.__structures.append_structure(atom)
            subatom = copy.deepcopy(atom)
            for etan in np.linspace(-0.05,0.05,21):

                for strains in range(len(subatom.strainList)):
                    subatom = copy.deepcopy(atom)
                    subatom.distort(eta=etan, strainType_index = strains)
                    self.__structures.append_structure(subatom)
        
        self.__structures.write_structures(self.__structures)
        print [struct[1].path for struct in self.__structures.get_structures().items()]
        #return self.__structures
    
    def generate_supercells(self):
        
        dirnames = [struct[1].path for struct in self.__structures.get_structures().items()]
        rootdir = os.getcwd()
        for d in dirnames:
            os.chdir(d)
            os.system('/home/MCL/t.dengg/bin/phonopy-1.7.4/bin/phonopy  -d --dim="5 5 5" -c POSCAR')
            os.system('cp POSCAR POSCAR-p')
            os.system('mv SPOSCAR POSCAR')
            os.system('pwd')
            os.chdir(rootdir)
            
    
    def distortions(self, volume, etamax, N):
        return self.__structures
    
class Postprocess(ECs):
    def __init__(self):
        super(Postprocess, self).__init__()
        self.__allstructures = None
        
    def set_energies(self):
        ec = ECs('vasp',True)
        ec.set_structures()
        ec.set_gsenergy()
        #ec.set_fenergy()
        self.__allstructures = ec.get_structures()
        
        
    def eq_volumes(self, T_array):
        
        return self.__vmin, T_array
    
    def rewrite_structures(self):
        structures = []
        self.__V=[]
        for (a,b,v) in self.__allstructures.keys():
            if not v in self.__V:
                self.__V.append(v)
        self.__V = sorted(self.__V)
        i=0
        
        for v in self.__V:
            
            structures.append({})
            for (a,b,c) in self.__allstructures.keys():
                if c==v:
                    
                    structures[i][(a,b)] = self.__allstructures[(a,b,v)]
            
            f=open('./%s/structures.pkl'%v,'w')
            pickle.dump(structures[i],f)
            f.close()
            i+=1 
            
        self.__structures = structures
        return 
    
    def ECs(self):
        ec = ECs('vasp',True)
        i=0
        for structure in self.__structures:
            ec = ECs('vasp',True)
            os.chdir(str(self.__V[i]))
            # T-dependent EC's ############
            T = structure.values()[0].T
            ec.T = T
            ec.set_structures()
            #print ec.get_structures()
            """
            eta=[]
            scale=[]
            gsenergy=[]

            for atom in ec.get_atomsByStraintype('01'):
                eta.append(atom.eta)
                scale.append(atom.scale)
                gsenergy.append(atom.gsenergy)
            print eta,gsenergy
            a1,b1,c1,d1 = ANALYTICS(eta,gsenergy).phist(50, 100, 0.0000001)
            eta=[]
            scale=[]
            gsenergy=[]
            for atom in ec.get_atomsByStraintype('08'):
                eta.append(atom.eta)
                scale.append(atom.scale)
                gsenergy.append(atom.gsenergy)
            print eta,gsenergy
            a2,b2,c2,d2 = ANALYTICS(eta,gsenergy).phist(50, 100, 0.0000001)
            eta=[]
            scale=[]
            gsenergy=[]
            for atom in ec.get_atomsByStraintype('23'):
                eta.append(atom.eta)
                scale.append(atom.scale)
                gsenergy.append(atom.gsenergy)
            print eta,gsenergy
            a3,b3,c3,d3 = ANALYTICS(eta,gsenergy).phist(50, 100, 0.0000001)
            
            ec.fitorder=[a1[1],a2[1],a3[1]]
            ec.etacalc=[str(a1[0]),str(a2[0]),str(a3[0])]
            """
            ec.fitorder = [6,6,6]
            ec.etacalc = ['0.05','0.05','0.04']
            ec.set_analytics()
            ecs = ec.get_ec()
            cvs = ec.get_CVS()
            
            #ec.plot_cvs()
            #ec.plot_cvs('E0')
            #ec.plot_cvs('Fvib')
            
            #ec.plot_2nd()
            #ec.plot_2nd('E0')
            #ec.plot_2nd('Fvib')
            
            #ec.plot_energy()
            #ec.plot_energy(mod='E0')
            #ec.plot_energy(mod='Fvib')
            
            f=open('ECs.pkl','w')
            pickle.dump(ecs,f)
            f.close()
            f=open('CVS.pkl','w')
            pickle.dump(cvs,f)
            f.close()
            ################################
            
            # T=0K EC's ####################
            ec = ECs('vasp')
            ec.set_structures()
            ec.fitorder=[6,6,6]
            ec.etacalc=['0.05','0.05','0.04']
            ec.set_analytics()
            ecs = ec.get_ec()
            cvs = ec.get_CVS()
            f=open('ECs_T0.pkl','w')
            pickle.dump(ecs,f)
            f.close()
            f=open('CVS_T0.pkl','w')
            pickle.dump(cvs,f)
            f.close()
            ################################
            
            
            os.chdir('../')
            i+=1
        return
    
    def interpolate(self):
        return self.__interpolated
    