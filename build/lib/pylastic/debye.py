import json
import numpy as np
from scipy.optimize import fmin, brent
from math import factorial as fc
class Debye():
    
    def __init__(self):
        self.__kb=8.617332478*10.**(-5.) #eV/K
        self.__h =4.135667516*10.**(-15.) #eV s
        self.__path = './'
        self.T = 1.
        self.__fitorder = 2
        self.m = 183.84 * 1.66053892173 * 10.**(-27.)#Atomic weight in kg
        self.__mod='X/B-fit'
        self.__E_fname = '/Energy_vol_W.dat'
        self.__C_fname = '/Cij_0.json'
    
    def get_mod(self):
        return self.__mod
    
    def set_mod(self, mod):
        print 'setting mod'
        self.__mod =  mod
    
    def get_fitorder(self):
        return self.__path
    
    def set_fitorder(self, fitorder):
        print 'setting fitorder'
        self.__path = fitorder
        
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
    
    def get_E0(self):
        return self.__E0
    def set_E0(self,E0):
        self.__E0 = E0
        
    def get_V(self):
        return self.__V
    def set_V(self,V):
        self.__V = V   
        
    def calculate_moduli(self, scale):
        Cij_V = self.__Cij_dic[scale]['SM']
        C = np.zeros((6,6))
        for j in range(6):
            for i in range(6):
                C[i,j] = Cij_V[i+6*j]
        BV = (C[0,0]+C[1,1]+C[2,2]+2*(C[0,1]+C[0,2]+C[1,2]))/9
        GV = ((C[0,0]+C[1,1]+C[2,2])-(C[0,1]+C[0,2]+C[1,2])+3*(C[3,3]+C[4,4]+C[5,5]))/15
        EV = (9*BV*GV)/(3*BV+GV)
        
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
        
        return D
    
    def debye_function_exp(self, x):
        Brillouin = [1./6.,-1./30.,1./42.,-1./30.,5./66.,-691./2730.,7./6.,-3617./510.,43867./798,-174611./330.]
        D2=0.
        for k in range(len(Brillouin)):
            D2 = D2+Brillouin[k]/( (2*(k+1)+3)*fc(2*(k+1)) )*x**(2*(k+1))
        return 1- 3./8.*x + 3*D2
    
    def debye_T(self, x):
        
        self.__EoB = []
        self.__GoB = []
        self.__G = []
        self.__E = []
        self.__B = []
        self.__V = []
        
        dic = self.get_Cij()
        for scale in sorted(dic.keys()):
            (a, b, c, d, e) = self.calculate_moduli(scale)
            self.__EoB.append(c)
            self.__GoB.append(b)
            self.__B.append(a)
            self.__G.append(d)
            self.__E.append(e)
            self.__V.append(float(scale)**3./2.*10.**(-30.))
        
        
        c1= np.polyfit(self.__V, self.__EoB, self.__fitorder)
        p_EoB = np.poly1d(c1)
           
        c2= np.polyfit(self.__V, self.__GoB, self.__fitorder)
        p_GoB = np.poly1d(c2)
        
        c3= np.polyfit(self.__V, self.__B, self.__fitorder)
        p_B = np.poly1d(c3)
        
        rho = self.m/x
        
        c4= np.polyfit(self.__V, self.__E, self.__fitorder)
        p_E = np.poly1d(c4)
        
        c5= np.polyfit(self.__V, self.__G, self.__fitorder)
        p_G = np.poly1d(c5)
        
        Const = self.__h/self.__kb* (3./(4.*np.pi))**(1./3.)
        if self.__mod=='X/B-fit':
            theta = Const * ( 1./3.*(p_EoB(x))**(-3./2.) + 2./3.*(p_GoB(x))**(-3./2.) )**(-1./3.) * (p_B(x)*10**(9.)/rho)**(1./2.) * x**(-1./3.)
        else:
            theta = Const * ( 1./3.*(p_E(x)/p_B(x))**(-3./2.) + 2./3.*(p_G(x)/p_B(x))**(-3./2.) )**(-1./3.) * (p_B(x)*10**(9.)/rho)**(1./2.) * x**(-1./3.)

        return theta
    
    def free_energy(self,x):
        E0 = self.get_gsenergy()
        
        T_Deb = self.debye_T(x)
        
        V0=[float(l)**3./2.*10.**(-30.) for l in self.__l]
        c1= np.polyfit(V0, E0, self.__fitorder)
        p_E0 = np.poly1d(c1)
        #print self.debye_function(self.debye_T(x)/self.T),self.debye_T(x),self.T
        return (p_E0(x) + ( self.debye_function(T_Deb/self.T) + 3.*np.log(1.-np.exp(-T_Deb/self.T)) ) * self.__kb * self.T - 9./8.*self.__kb*T_Deb)
    
    def free_energy_vib(self,x):
        E0 = self.get_gsenergy()
        
        
        #print self.__GoB
        #print -(  self.debye_function(self.debye_T(x)/self.T) + 3.*np.log(1.-np.exp(-self.debye_T(x)/self.T)) ) * self.__kb * self.T,9./8.*self.__kb*self.debye_T(x)
        return +(  self.debye_function(self.debye_T(x)/self.T) + 3.*np.log(1.-np.exp(-self.debye_T(x)/self.T)) ) * self.__kb * self.T + 9./8.*self.__kb*self.debye_T(x)
    
        
    
    def optimization(self):
        #return brent(self.free_energy, brack=(15.5*10.**(-30.),17.*10.**(-30.)))
        return fmin(self.free_energy, 15.6*10.**(-30.))
    
    def find_min(self, listx, listy):
        minval = min(listy)
        minindex = listy.index(minval)
        m = listx[minindex]
        return m
    
    E0 = property(fget=get_E0, fset=set_E0)
    V = property(fget=get_V, fset=set_V)  
    path = property(fget=get_path, fset=set_path)
    fitorder = property(fget=get_fitorder, fset=set_fitorder)
    E_fname = property(fget=get_E_fname, fset=set_E_fname)
    C_fname = property(fget=get_C_fname, fset=set_C_fname)
    mod = property(fget=get_mod, fset=set_mod)