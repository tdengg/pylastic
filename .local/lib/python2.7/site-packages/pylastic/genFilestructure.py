import os
import sys
import numpy as np


class Filestructure(object):
    
    def __init__(self):
        self.__path == None
        
    def generate(self):

        dic = {}
        n = 1
        for distype in self.__strainList:
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

self.__strainList
self.__eta
self.__NoP  