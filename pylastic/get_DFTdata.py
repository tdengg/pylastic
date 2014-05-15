import lxml.etree as et

class VASP(object):
    def __init__(self):
        self.__vfile = None
        self.__cellsize = 1
        self.__ERange = (-3000.,-1000.)
    ##    
    def set_outfile(self, vfile):
        self.__vfile = vfile
        self.__vasprun = et.parse(self.__vfile)
    
    def get_outfile(self):
        return self.__vasprun
    ##
    
    def set_gsEnergy(self):
        elem = self.__vasprun.xpath("//scstep[last()]/energy/i[@name = 'e_fr_energy']")
        self.__gsEnergy = float(elem[0].text)
        
    #####################################
    def set_cellsize(self, cells):
        self.__cellsize = cells
        
    def get_cellsize(self):
        return self.__cellsize
    
    def set_ERange(self, ERange):
        """Set expected energy range to filter VASP output: as tuple (emin,emax)"""
        self.__ERange = ERange
    def set_gsEnergy_DFPT(self):
        """Compensate for VASP's messy output"""
        elem = self.__vasprun.xpath("//scstep/energy/i[@name = 'e_fr_energy']")
        
        allengys = []
        for k in elem:
            try:
                allengys.append(float(k.text))
            except:
                allengys.append(0.)
        
        trueengys = []
        for engy in allengys:
            if engy < self.__ERange[1] and engy > self.__ERange[0]: trueengys.append(engy)
        self.__gsEnergy = trueengys[-1]/self.__cellsize   
    ####################################
    
    
    def get_gsEnergy(self):
        return self.__gsEnergy
    
    ##
    def set_freeEnergy(self, T):
        """Process PHONOPY free energy output"""
        g = open(self.vfile)
        self.__freeEnergy = float(g.readlines()[T].split()[1])/96.47244
        g.close()
        
    def get_freeEnergy(self):
        return self.__freeEnergy
    ##