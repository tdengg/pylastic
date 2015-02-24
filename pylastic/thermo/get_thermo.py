

class Phonon(object):
    def __init__(self):
        self.__dir = '.'
        self.__T = None
        self.__F = None
    
    def set_path(self, path):
        self.__dir = path.rstrip('/')
    
    def get_path(self):
        return self.__dir
    
    def set_T(self, T):
        if not self.__T: self.set_F()
        if not type(self.__T) == 'list': self.set_F()
        if T in self.__T:
            self.__F = self.__F[self.__T.index(T)]
            self.__T = T
        else: print "Warning: No calculated values for given temperature - setting it to %f K"%self.__T
            
    def get_T(self):
        return self.__T
    
    def set_F(self):
        f = open(self.__dir + '/F_TV', 'r')
        fdata = f.readlines()
        f.close()
        T = []
        F = []
        for i in fdata:
            T.append(float(i.split()[0]))
            F.append(float(i.split()[1])/96.47244)
            
        self.__F = F
        self.__T = T
        
    def get_F(self):
        return self.__T, self.__F
    
    def set_phononDOS(self):
        f = open(self.__dir + '/total_dos.dat', 'r')
        dos = f.readlines()
        f.close()
        
        f = []
        di = []
        i=0
        for line in dos:
            if i == 0: 
                i+=1
                continue
            
            di.append(float(line.split()[1]))
            f.append(float(line.split()[0]))
            i+=1
        self.__dos = di
        self.__f = f
        
    def get_phononDOS(self):
        return self.__f, self.__dos
    
    