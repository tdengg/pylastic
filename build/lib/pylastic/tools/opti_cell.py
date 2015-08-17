import numpy as np
from pylastic.elatoms import ElAtoms, Structures
from pylastic.postprocess import ECs

class HCP(object):
    def __init__(self, cod='vasp'):
        self.__cod=cod
        if self.__cod=='vasp': 
            from pylastic.io.vasp import POS
            
        
        poscar = POS('POSCAR').read_pos()
        
        structures = Structures(self.__cod)

        ## Generate distorted structures and add them to structures object: ##
        atom = ElAtoms(self.__cod)
        atom.poscarToAtoms(poscar)
        
        
        i=0
        imax = 100
        while i<imax:
            
            atom = ElAtoms(self.__cod)
            atom.poscarToAtoms(poscar)
            for scale in np.linspace(amin,amax,N):

                atom = ElAtoms(self.__cod)
                atom.poscarToAtoms(poscar)
                atom.scale = scale
                structures.append_structure(atom)
            
            ec = ECs(self.__cod)
            ec.set_gsenergy()
            
            eta=[]
            scale=[]
            gsenergy=[]
            for atom in ec.get_atomsByStraintype(None):
                eta.append(atom.eta)
                scale.append(atom.scale)
                gsenergy.append(atom.gsenergy)
            
            
            #polyfit and set new min
            aopti = min
            for coa in range():
                #set structure
                
                #calculate
            #polyfit and set new
            
            
            delta= abs(currcoamin-coamin)
            coamin = min 
            
            
            i+=1
        return
    
    def 
        