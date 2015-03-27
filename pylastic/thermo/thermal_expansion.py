import numpy as np
from pylastic.io.espresso import POS
import copy
import os

class Setup(object):
    def __init__(self):
        self.__N = None
        self.__delta = None
        self.__fname = 'si.scf.in'
        self.__code = 'espresso'
        
    def eos(self, N, delta):
        if self.__code == 'espresso':
            P = POS(self.__fname)
            dic = P.read_in()
            
            for i in range(N):
                dicnew = copy.deepcopy(dic)
                dicnew['vlatt_1'] = np.array(dicnew['vlatt_1'])+np.array(dicnew['vlatt_1'])*delta*(i-N/2.)
                dicnew['vlatt_2'] = np.array(dicnew['vlatt_2'])+np.array(dicnew['vlatt_2'])*delta*(i-N/2.)
                dicnew['vlatt_3'] = np.array(dicnew['vlatt_3'])+np.array(dicnew['vlatt_3'])*delta*(i-N/2.)
                try:
                    os.mkdir(str(delta*(i-N/2.)))
                except:
                    print 'Directory %s exists'%str(delta*(i-N/2.))
                P.write_in(dicnew,'./%s/espresso.in'%str(delta*(i-N/2.)))
            print dic
        return