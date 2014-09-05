""" ElaStic standard setup
"""

import numpy as np
import os
import json

import distort
import spacegroup
import io.vasp

class VASP(object):
    def __init__(self, eta=0.05, method='Energy', order=2, NoP=11, fname='POSCAR'):
        self.__eta = eta
        self.__method = method
        self.__order = order
        self.__NoP = NoP
        self.__fname = fname
        self.__paths = []
    
    def gen_filestruct(self):
        o_poscar = io.vasp.POS(self.__fname)
        poscar = o_poscar.read_pos()                 # Generate sgroup.out from POSCAR file                    ###MOD
        
        
        ###### Structure ####################
        bv = [poscar["vlatt_1"],poscar["vlatt_2"],poscar["vlatt_3"]]
        M_old= np.array(bv)
        D    = np.linalg.det(M_old)
        V0   = abs(poscar['scale']**3*D)
        #####################################
        
        sg = spacegroup.Sgroup(poscar)                      #spacegroup instance
        dist = distort.Distort()                            #distort instance
        
        dist.sgn = sg.sgn                                   #get and set spacegroup number
        
        dist.set_strainList()                               # set strain list according to space group number
        dist.set_V0(V0)
        dic = {}
        n = 1
        for distype in dist.strainList:
            dic[str(n)] = {}
            dst=1
            os.mkdir(os.getcwd()+'/Dst%.2d'%n)
            for eta in np.linspace(-self.__eta,self.__eta, self.__NoP):
                dist.eta = eta
                dist.set_strainType(distype)
                dist.set_defMatrix()
                def_mat = dist.get_defMatrix()
                M_new = np.dot(M_old, def_mat)
                
                
                pname = os.getcwd()+"/Dst%(type).2d/Dst%(type).2d_%(numb).2d/"%{'type':int(n),'numb':int(dst)}
                for j in range(3):           
                    poscar["vlatt_%s"%(j+1)] = [M_new[j,0],M_new[j,1],M_new[j,2]]
                os.mkdir(pname)
                o_poscar.write_pos(poscar, pname+'POSCAR')
                os.system('cp INCAR '+pname)
                os.system('cp POTCAR '+pname)
                os.system('cp KPOINTS '+pname)
                info = dist.get_info()
                self.__paths.append(pname)
                dic[str(n)][str(dst)] = info
                
                
                dst+=1
            n+=1
        
        out = json.dumps(dic, indent=4)
        f=open('info.json','w')
        f.write(out)
        f.close()
        return
    
    def calc(self):
        root = os.getcwd()
        if not self.__executable:
            raise ValueError('Please define location of vasp executable!')
        for path in self.__paths:
            os.chdir(path)
            print 'starting vasprun for ',path
            os.system(self.__executable)
            os.chdir(root)
        return
    
    def set_executable(self, executable):
        self.__executable = executable
        

if __name__ == '__main__':
    calc = VASP(0.05,'Energy',2,11,'POSCAR')
    calc.gen_filestruct()
    calc.set_executable('/home/t.dengg/bin/vasp/vasp.5.3/vasp')
    calc.calc()