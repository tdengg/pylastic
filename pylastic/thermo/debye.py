import json
import numpy as np
from pylastic.tools.eos import Birch
from scipy.optimize import fmin, brent
from math import factorial as fc

import matplotlib.pyplot as plt

class Debye(object):
    
    def __init__(self):
        self.__kb=8.6173303*10.**(-5.) #eV/K
        self.__h =4.135667662*10.**(-15.) #eV s
        self.__path = './'
        self.T = 1.
        self.__fitorder_EC = 3
        self.__fitorder_EOS = 6
        self.m = 183.84 * 1.66053892173 * 10.**(-27.)#Atomic weight in kg
        self.__mod='X/B-fit'
        self.__E_fname = '/Energy_vol_W.dat'
        self.__C_fname = '/Cij_0.json'
        self.__artificial_deformation=np.zeros((6,6))
        self.__artdef=False
        self.__elastic = True
        self.__V0 = None
        self.__lt = False
        self.__natom=2.
        #self.__fout_EC = open('out_EC','a')
    
    def set_artificial_deformation(self, eta=0.1 ,typeis='tetragonal'):
        self.__artdef = True
        if typeis=='tetragonal':
            print 'Setting artificial deformation to tetragonal'
            
            self.__artificial_deformation[0,0] = eta
            self.__artificial_deformation[1,1] = eta
            self.__artificial_deformation[2,2] = eta
            #self.__artificial_deformation[0,1] = eta
            #self.__artificial_deformation[1,0] = eta
            #self.__artificial_deformation[2,0] = eta
            #self.__artificial_deformation[3,3] = eta
            #self.__artificial_deformation[4,4] = eta
            #self.__artificial_deformation[5,5] = eta
    
    def get_numatom(self):
        return self.__natom
    
    def set_numatom(self, natom):
        self.__natom = natom
    
    def get_mod(self):
        return self.__mod
    
    def set_mod(self, mod):
        print 'setting mod: %s'%mod
        self.__mod =  mod
        
    def get_lt(self):
        return self.__lt
    
    def set_lt(self, lt):
        print 'setting low temperature correction: %s'%lt
        self.__lt =  lt
        
    def get_elastic(self):
        return self.__elastic
    
    def set_elastic(self, elastic):
        print 'elastic = %s'%elastic
        self.__elastic =  elastic 
    
    def get_fitorder_EC(self):
        return self.__fitorder_EC
    
    def set_fitorder_EC(self, fitorder):
        print 'setting fitorder'
        self.__fitorder_EC = fitorder
        
    def get_fitorder_EOS(self):
        return self.__fitorder_EOS
    
    def set_fitorder_EOS(self, fitorder):
        print 'setting fitorder'
        self.__fitorder_EOS = fitorder        
        
    def get_E_fname(self):
        return self.__E_fname
    
    def set_E_fname(self, E_fname):
        print 'setting E_fname'
        self.__E_fname =  E_fname
        
    def get_C_fname(self):
        return self.__C_fname
    
    def set_C_fname(self, C_fname):
        print 'setting C_fname'
        self.__C_fname =  C_fname
    
    def get_path(self):
        return self.__path
    
    def set_path(self, path):
        print 'setting path'
        self.__path = path
    
    def get_gsenergy(self):
        f=open(self.__path+self.__E_fname)
        lines = f.readlines()
        f.close()
        self.__E0 = [float(line.split()[0]) for line in lines]
        self.__l = [float(line.split()[1]) for line in lines]
        return self.__E0
    
    def get_Cij(self):
        f=open(self.__path + self.__C_fname)
        self.__Cij_dic = json.load(f)
        f.close()
        return self.__Cij_dic
    
    def calc_B(self):
        birch = Birch(self.__V,self.__E0).fit()
        B = birch.out1
        print B
        return B
    
    def get_E0(self):
        return self.__E0
    def set_E0(self,E0):
        self.__E0 = E0
        
    def get_V(self):
        return self.__V
    def set_V(self,V):
        self.__V = V   
    
    def get_V0(self):
        return self.__V0
    def set_V0(self,V0):
        self.__V0 = V0
    
    def calculate_moduli(self, scale):
        Cij_V = np.array(self.__Cij_dic[scale]['SM'])
        C = np.zeros((6,6))
        if not Cij_V.shape==(6,6):
            for j in range(6):
                for i in range(6):
                    C[i,j] = Cij_V[i+6*j]
        else:
            C=Cij_V
        if self.__artdef:   
            #print 'Modifying elastic tensor.'
            C = C+self.__artificial_deformation
        BV = (C[0,0]+C[1,1]+C[2,2]+2.*(C[0,1]+C[0,2]+C[1,2]))/9.*1.4
        GV = ((C[0,0]+C[1,1]+C[2,2])-(C[0,1]+C[0,2]+C[1,2])+3.*(C[3,3]+C[4,4]+C[5,5]))/15.
        #EV = (9.*BV*GV)/(3.*BV+GV)
        #### calculate p-wave modulus ####
        EV = BV + 4./3.*GV
        #if self.__outmod == 'massiveoutput':
        #   self.__fout_EC.write(str(scale) +' '+ str(BV) +' '+ str(GV) +' '+ str(EV) + '\n')
        return BV, GV/BV, EV/BV, GV, EV
    
    def debye_function(self, x):
        exx = np.exp(-x)
        if x == 0: D=1.
        elif 0 < x <= 0.1:
            D = 1 - 0.375*x + x**2.*(0.05-5.952380953*10**(-4.)*x**2.)
        elif 0.1 < x <= 7.25:
            D = ( (((0.0946173*x-4.432582)*x+85.07724)*x-800.6087)*x+3953.632 ) / ( (((x+15.121491)*x+143.155337)*x+682.0012)*x+3953.632 )
        elif x > 7.25:
            N=int(25./x)
            
            D=0.
            D2=1.
            
            
            if N>1.:
                for i in range(1,N+1):
                    
                    DS = i
                    
                    x3 = DS*x
                    D2 = D2*exx
                    D = D + D2*( 6. + x3*(6. + x3*(3.+x3)) )/DS**4.
                
                    
            
                
            D = 3.*(6.493939402-D)/(x**3.)
        else: 
            print 'ERROR: Debye function out of bounds: D(%s)!'%(x)
        return D
    
    def debye_function_exp(self, x):
        Brillouin = [1./6.,-1./30.,1./42.,-1./30.,5./66.,-691./2730.,7./6.,-3617./510.,43867./798,-174611./330.]
        D2=0.
        for k in range(len(Brillouin)):
            D2 = D2+Brillouin[k]/( (2*(k+1)+3)*fc(2*(k+1)) )*x**(2*(k+1))
        return 1.- 3./8.*x + 3.*D2
    
    def debye_T(self, x):
        
        self.__EoB = []
        self.__GoB = []
        self.__G = []
        self.__E = []
        self.__B = []
        self.__V = []
        if self.__elastic:
            dic = self.get_Cij()
            for scale in sorted(dic.keys()):
            
                (a, b, c, d, e) = self.calculate_moduli(scale)
                self.__EoB.append(c)
                self.__GoB.append(b)
                self.__B.append(a)
                self.__G.append(d)
                self.__E.append(e)
                self.__V.append(float(scale)**3./self.__natom*10.**(-30.))
            
            
            c1= np.polyfit(self.__V, self.__EoB, self.__fitorder_EC)
            p_EoB = np.poly1d(c1)
               
            c2= np.polyfit(self.__V, self.__GoB, self.__fitorder_EC)
            p_GoB = np.poly1d(c2)
            
            c3= np.polyfit(self.__V, self.__B, self.__fitorder_EC)
            p_B = np.poly1d(c3)
            
            rho = self.m/x
            
            c4= np.polyfit(self.__V, self.__E, self.__fitorder_EC)
            p_E = np.poly1d(c4)
            
            c5= np.polyfit(self.__V, self.__G, self.__fitorder_EC)
            p_G = np.poly1d(c5)
        
        else:
            self.__mod = 'prefactor'
            for scale in sorted(dic.keys()):
                self.calc_B()
            
            c3= np.polyfit(self.__V, self.__B, self.__fitorder_EC)
            p_B = np.poly1d(c3)
            ############# Lorenz ############
            
            a= 17.343  
            b= -363.343
            c= 881.843 
            d= 1077.8 
            #p_B = lambda x: (10.0/9*b*x**(-8.0/3)+28.0/9*c*x**(-10.0/3)+6*d*x**(-4))*x*160.22
            #################################
        if self.__lt: 
            lowT_correction = (x/self.__V0)**(1./3.)
            
        else: lowT_correction = 1.
        
        Const = self.__h/self.__kb* (3./(4.*np.pi))**(1./3.)
        if self.__mod=='X/B-fit':
            print ( 1./3.*(p_EoB(x))**(-3./2.) + 2./3.*(p_GoB(x))**(-3./2.) )**(-1./3.)
            theta = Const * ( 1./3.*(p_EoB(x))**(-3./2.) + 2./3.*(p_GoB(x))**(-3./2.) )**(-1./3.) * (p_B(x)*10**(9.)/rho)**(1./2.) * x**(-1./3.)*lowT_correction
        elif self.__mod=='prefactor': 
            theta = Const * 0.617 * (p_B(x)*10**(9.)/rho)**(1./2.) * x**(-1./3.)*lowT_correction
        elif self.__mod=='debug':
            B0=321.924303414506
            V0=15.7103029969695*10.**(-30)
            rho0=self.m/V0
            gamma= 1.62296663327697
            theta0 = Const * 0.617 * (B0*10**9./rho0)**(1./2.) * V0**(-1./3.)*lowT_correction
            #print theta0,rho0,rho,V0,x,p_B(x)
            theta = theta0*(V0/x)**gamma
        else:
            theta = Const * ( 1./3.*(p_E(x)/p_B(x))**(-3./2.) + 2./3.*(p_G(x)/p_B(x))**(-3./2.) )**(-1./3.) * (p_B(x)*10**(9.)/rho)**(1./2.) * x**(-1./3.)*lowT_correction

        return theta
    
    def free_energy(self,x):
        self.__E0 = self.get_gsenergy()
        
        V0=[float(l)**3./self.__natom*10**(-30) for l in self.__l]
        c1= np.polyfit(V0, self.__E0, self.__fitorder_EOS)
        p_E0 = np.poly1d(c1)
        dp_E0 = np.polyder(p_E0)
        roots = np.roots(dp_E0)
        isreal = np.isreal(roots)
        i=0
        for val in isreal:
            if val: index = i
            
            i+=1
        
        self.__V0 = roots[index]
        
        self.T_Deb = self.debye_T(x)
        ##########Lorenz############
        a= 17.343  
        b= -363.343
        c= 881.843 
        d= 1077.8 
        #p_E0 = lambda X: a+b*(X*10**(30))**(-2.0/3)+c*(X*10**(30))**(-4.0/3)+d*(X*10.**(30))**(-2.)
        
        #Vtest1=[15.*10**(-30),16.*10**(-30),17.*10**(-30)]
        Vtest=[15.,16.,17.]
        #plt.plot(Vtest,[p_E01(v) for v in Vtest])
        #plt.plot(Vtest,[p_E0(v) for v in Vtest1])

        #plt.show()
        ############################
        #print self.debye_function(self.debye_T(x)/self.T),self.debye_T(x),self.T
        #return (p_E0(x) - ( self.debye_function(self.T_Deb/self.T) - 3.*np.log(1.-np.exp(-self.T_Deb/self.T)) ) * self.__kb * self.T + 9./8.*self.__kb*self.T_Deb)
        return (p_E0(x) + (( -self.debye_function(self.T_Deb/self.T) + 3.*np.log(1.-np.exp(-self.T_Deb/self.T)) ) * self.__kb * self.T - 9./8.*self.__kb*self.T_Deb))

    def free_energy_vib(self,x):
        #E0 = self.get_gsenergy()
        
        
        #print self.__GoB
        #print -(  self.debye_function(self.debye_T(x)/self.T) + 3.*np.log(1.-np.exp(-self.debye_T(x)/self.T)) ) * self.__kb * self.T,9./8.*self.__kb*self.debye_T(x)
        return +(   self.debye_function(self.debye_T(x)/self.T) + 3.*np.log(1.-np.exp(-self.debye_T(x)/self.T)) ) * self.__kb * self.T + 9./8.*self.__kb*self.debye_T(x)
    
        
    def optimization(self):
        #return brent(self.free_energy, brack=(15.5*10.**(-30.),17.*10.**(-30.)))
        return fmin(self.free_energy, 15.6)
    
    def find_min(self, listx, listy):
        minval = min(listy)
        minindex = listy.index(minval)
        m = listx[minindex]
        return m
    
    E0 = property(fget=get_E0, fset=set_E0)
    V = property(fget=get_V, fset=set_V)  
    V0 = property(fget=get_V0, fset=set_V0)  
    path = property(fget=get_path, fset=set_path)
    fitorder_EC = property(fget=get_fitorder_EC, fset=set_fitorder_EC)
    fitorder_EOS = property(fget=get_fitorder_EOS, fset=set_fitorder_EOS)
    E_fname = property(fget=get_E_fname, fset=set_E_fname)
    C_fname = property(fget=get_C_fname, fset=set_C_fname)
    mod = property(fget=get_mod, fset=set_mod)
    lt = property(fget=get_lt, fset=set_lt)
    elastic = property(fget=get_elastic, fset=set_elastic)
    natom = property(fget=get_numatom, fset=set_numatom)