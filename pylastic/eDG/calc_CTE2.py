import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.optimize import fmin
from scipy.integrate import quad
from scipy.optimize import brentq

from pylastic.tools.eos import Birch
from pylastic.tools.birch import EOS
import os


kb=np.float64(1.3866*10.**(-23.))
N=6.022*10**23
h=np.float64(6.626*10.**(-34.)/(2.*np.pi))

eV=6.241509126*10**18.#1.602*10**(-19.)

A=10000.

ax1=plt.subplot(131)
ax2=plt.subplot(132)
ax3=plt.subplot(133)

print os.getcwd()

def debye_function(x):
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

def alpha(trange,minl):
    alph = []
    
    for i in range(len(trange)-1):
        
        Alpha = (minl[i+1]-minl[i])/(trange[i+1]-trange[i])/minl[0]*10.**6.
        if Alpha >= 0:
            alph.append(Alpha/3.)
        else: 
            alph.append(0)
    trange.pop(0)
    
    return trange,alph


def Integrand(x,h,kb,c,T, V,V0,B,A):
    
    Vs=V
    #return kb*T*Vs**(1./3.)*A**2./(2.*np.pi**2.*c)*(np.exp(x*(Vs)**(1./3.)/(c*A))-1.)**2.*np.exp(x*(Vs)**(1./3.)/(c*A))*( h*x/(kb*T)*(1./(np.exp(h*x/(kb*T))-1)+1.)-np.log(np.exp(h*x/(kb*T))-1.) )
    #return Vs**(1./3.)*A**2./(2.*np.pi**2.*c)*h*x*(np.exp(x*(Vs)**(1./3.)/(c*A))-1.)**2.*np.exp(x*(Vs)**(1./3.)/(c*A))*(1./(np.exp(h*x/(kb*T))-1)+1./2.)
    return -kb*T*Vs**(1./3.)*A**2./(2.*np.pi**2.*c)*(np.exp(x*(Vs)**(1./3.)/(c*A))-1.)**2.*np.exp(x*(Vs)**(1./3.)/(c*A))*np.log(1.-np.exp(-h*x/(kb*T)))#-B*(V0-V)
    
def IntegrandD(x,h,kb,V,c,T,V0,B):
    Vs=V
    return -Vs*kb*T/(2.*np.pi**2.*c**3.)*x**2.*np.log(1.-np.exp(-h*x/(kb*T)))#-B*(V-V0)
    #return kb*T*x**2.*Vs/(2.*np.pi**2.*c**3.)*( h*x/(kb*T)*(1./(np.exp(h*x/(kb*T))-1)+1.)-np.log(np.exp(h*x/(kb*T))-1.) )


def read_BGL(fname):
    f=open(fname)
    Cij_dic = json.load(f)
    f.close()
    
    a=[]
    E=[]
    B=[]
    G=[]
    G2=[]
    Cp=[]
    C11=[]
    C12=[]
    C44=[]
    V0_c=[]
    A=[]
    for scale in sorted(Cij_dic.keys()):
        Cij_V = np.array(Cij_dic[scale]['SM'])
        C = np.zeros((6,6))
        if not Cij_V.shape==(6,6):
            for j in range(6):
                for i in range(6):
                    C[i,j] = Cij_V[i+6*j]
        else:
            C=Cij_V

        BV = (C[0,0]+C[1,1]+C[2,2]+2.*(C[0,1]+C[0,2]+C[1,2]))/9.
        GV = ((C[0,0]+C[1,1]+C[2,2])-(C[0,1]+C[0,2]+C[1,2])+3.*(C[3,3]+C[4,4]+C[5,5]))/15.
        GV2 = ((C[0,0]+C[1,1]+C[2,2])-(C[0,1]+C[0,2]+C[1,2])+3.*(C[3,3]+C[4,4]+C[5,5]))/15.
        
        #### calculate p-wave modulus ####
        EV = BV + 4./3.*GV

    
        B.append(BV)
        E.append(EV)
        G.append(GV)
        G2.append(GV2*10.**9.)
        Cp.append((C[0,0]-C[0,1])/2.*10.**9.)
        C11.append((C[0,0])*10.**9.)
        C12.append((C[0,1])*10.**9.)
        C44.append((C[3,3])*10.**9.)
        V0_c.append(np.float64((float(scale)**3./2.)*10**(-30.)))
        a.append(float(scale))
        A.append((C[3,3])/((C[0,0]-C[0,1])/2.))


    
    
    
    coeff = np.polyfit(a,B,2)
    p_B = np.poly1d(coeff)
    dp_B = np.polyder(p_B)
    
    coeff = np.polyfit(a,G,2)
    p_G = np.poly1d(coeff)
    dp_G = np.polyder(p_G)
    
    
    coeff = np.polyfit(a,E,2)
    p_E = np.poly1d(coeff)
    dp_E = np.polyder(p_E)
    
    coeff = np.polyfit(a,Cp,2)
    p_Cp = np.poly1d(coeff)
    dp_Cp = np.polyder(p_Cp)
    return B, G, E, p_B,p_G,p_E, dp_B, dp_G, dp_E

def get_emin(fname):
    f = open(fname)
    lines = f.readlines()
    V0 = [(float(l.split()[0])*10.**(-30.)) for l in lines]
    E0 = [float(l.split()[1]) for l in lines]
    f.close()
    eos = EOS(np.array([np.float64(v*10.**30.) for v in V0]),E0)
    
    eos.fit_birch()
    par = eos.par
    a0=(2.*par[0])**(1./3.)
    v0=par[0]*10.**(-30.)
    return a0, v0
    
###########################################################################################

B_c=[]
dB_c=[]
G_c=[]
dG_c=[]
E_c=[]
dE_c=[]

Vmin_c=[]
a_c=[]
rho_c=[]

kb=np.float64(1.3866*10.**(-23.))
N=6.022*10**23
h=np.float64(6.626*10.**(-34.)/(2.*np.pi))

conz=['00','03','06','09','12','18','25','50']
m_c = [183.85 * 1.66053892173 * 10.**(-27.),183.92 * 1.66053892173 * 10.**(-27.),183.99 * 1.66053892173 * 10.**(-27.),184.06 * 1.66053892173 * 10.**(-27.),184.13 * 1.66053892173 * 10.**(-27.),184.27 * 1.66053892173 * 10.**(-27.),184.43 * 1.66053892173 * 10.**(-27.),185.02 * 1.66053892173 * 10.**(-27.)]


i=0
for (c) in (conz):
    
    a0, v0 = get_emin('eos%s.dat'%c)
    
    B, G, E,p_B,p_G,p_E, dp_B, dp_G, dp_E = read_BGL('Cij_%s.json'%c)
    
    
    B_c.append(p_B(a0))
    dB_c.append(dp_B(a0))
    
    G_c.append(p_G(a0))
    dG_c.append(dp_G(a0))
    
    E_c.append(p_E(a0))
    dE_c.append(dp_E(a0))
    
    Vmin_c.append(v0)
    a_c.append(a0)
    rho_c.append(v0/m_c[i])
    
    i+=1


coeff = np.polyfit(map(float,conz),dB_c,1)
p_dBc = np.poly1d(coeff)
coeff = np.polyfit(map(float,conz),dG_c,1)
p_dGc = np.poly1d(coeff)
coeff = np.polyfit(map(float,conz),dE_c,1)
p_dEc = np.poly1d(coeff)

ax2.plot(map(float,conz), B_c)
ax2.plot(map(float,conz), E_c)
ax2.plot(map(float,conz), G_c)

ax2.plot(map(float,conz), rho_c)

ax3.plot(map(float,conz),p_dBc(map(float,conz)))
ax3.plot(map(float,conz),p_dEc(map(float,conz)))
ax3.plot(map(float,conz),p_dGc(map(float,conz)))
ax3.plot(map(float,conz), dB_c,'bo')
ax3.plot(map(float,conz), dE_c,'go')
ax3.plot(map(float,conz), dG_c,'ro')
###########################################################################################
print a_c

xp = np.linspace(0,10000000000000,1000)


TEMP_plot = [100.,1000.,1400.,1800.,2200.]    
ALPHA1 = []
ALPHA2 = []
ALPHA3 = []
ALPHA4 = []
ALPHA5 = []
ALPHA1e = []
ALPHA2e = []
ALPHA3e = []
ALPHA4e = []
ALPHA5e = []
ALPHA1a = []
ALPHA2a = []
ALPHA3a = []
ALPHA4a = []
ALPHA5a = []
############################ W ########################################################
masses=[np.float64(183.85 * 1.66053892173 * 10.**(-27.)),np.float64(183.92 * 1.66053892173 * 10.**(-27.)),np.float64(183.99 * 1.66053892173 * 10.**(-27.)),np.float64(184.06 * 1.66053892173 * 10.**(-27.)),np.float64(184.13 * 1.66053892173 * 10.**(-27.)),np.float64(184.27 * 1.66053892173 * 10.**(-27.)),np.float64(184.43 * 1.66053892173 * 10.**(-27.)),np.float64(185.02 * 1.66053892173 * 10.**(-27.))]
concentrations=['00','03','06','09','12','18','25','50']
AAs=[160.,140.,120.,100.,80.,60.,40.,20.]
AAs=[1000.]*8
delta_theta=np.linspace(0.9,1.,8)
delta_theta=[1.]*8

for (AA,con,m) in zip(AAs,concentrations,masses):
    
    k=0
    f = open('eos%s.dat'%con)
    lines = f.readlines()
    V0 = [(float(l.split()[0])*10.**(-30.)) for l in lines]
    E0 = [float(l.split()[1]) for l in lines]
    f.close()
    eos = EOS(np.array([np.float64(v*10.**30.) for v in V0]),E0)
    eos.fit_birch()
    par = eos.par

    
    Bc=B_c
    Gc=G_c
    Ec=E_c
    dB=p_dBc(float(con))
    dG=p_dGc(float(con))
    dE=p_dEc(float(con))
    

    if k==0:
        f_out=open('a_T.dat','w')
        f_out.write('V_min at T=0K: %s \n\n'%Vmin_c[k])

    
    v0=Vmin_c[k]
    
    Lmin=[]
    LminD=[]
    Lminstd=[]
    Vexp=V0
    TEMP = np.linspace(1.,np.float64(2200),220)
    j=0
    for T in TEMP:
        V=[]
        f=[]
        fd=[]
        fstd=[]
        F=[]
        FD=[]
        Fstd=[]
        F_v=[]
        FD_v=[]
        Fstd_v=[]


        i=0
        for v in Vexp:
            Bv0 = Bc[k]
            Bv=Bc[k]+((v*2.)**(1./3.)-(v0*2.)**(1./3.))*10.**(10.)*dB
            Gv0 = Gc[k]
            Gv=Gc[k]+((v*2.)**(1./3.)-(v0*2.)**(1./3.))*10.**(10.)*dG
            Ev0 = Ec[k]
            Ev=Ec[k]+((v*2.)**(1./3.)-(v0*2.)**(1./3.))*10.**(10.)*dE
            #print Bv0,Bv
            rho = m/v
            rho0 = m/v0
            V.append(v)
            vs=v
            lowt=1.#(v/v0)**(1./3.)
            
            C=( 1./3.*((Ev0)/Bv0)**(-3./2.) + 2./3.*(Gv0/Bv0)**(-3./2.) )**(-1./3.)


            T_Deb = 2.*np.pi*h/kb* (3./(4.*np.pi))**(1./3.)* C * np.sqrt(Bv*10.**9./rho) * v**(-1./3.)*lowt

            A=AA#*(T_Deb/T)
            if j==0: c=np.float64(np.sqrt(Ev0*10.**9./rho0))
            else: 
                c=np.float64(np.sqrt(Ev*10.**9./rho))
                #c=C*np.sqrt(Bv*10.**9./rho)
                #c=C*np.sqrt(Ev0*Bv/Bv0*10.**9./rho)
                
            y = lambda x: (A**3./(2.*np.pi**2.)*( 1./3.*np.exp(3.*x*(vs)**(1./3.)/(c*A)) - np.exp(2.*x*(vs)**(1./3.)/(c*A)) + np.exp(x*(vs)**(1./3.)/(c*A)) - 1./3. ))-1.
            omega = brentq(y, 10**10., 10**15.)*delta_theta[k]
            theta = omega*h/kb
            

            ft1=-quad(Integrand,0,omega, args=(h,kb,c,T,vs,v0,Bv*10.**9.,A))[0]*eV
            fD1=-quad(IntegrandD,0,(6.*np.pi**2./vs)**(1./3.)*c, args=(h,kb,v,c,T,v0,Bv*10.**9.))[0]*eV

            if j==0: c=np.float64(np.sqrt(Gv0*10.**9./rho0))
            else: 
                c=np.float64(np.sqrt(Gv*10.**9./rho))
                #c=C*np.sqrt(Bv*10.**9./rho)
                #c=C*np.sqrt(Gv0*Bv/Bv0*10.**9./rho)
                
            y = lambda x: (A**3./(2.*np.pi**2.)*( 1./3.*np.exp(3.*x*(vs)**(1./3.)/(c*A)) - np.exp(2.*x*(vs)**(1./3.)/(c*A)) + np.exp(x*(vs)**(1./3.)/(c*A)) - 1./3. ))-1.
            omega = brentq(y, 10**10., 10**15.)*delta_theta[k]
            theta = omega*h/kb
            
            
            ft2=-quad(Integrand,0,omega, args=(h,kb,c,T,vs,v0,Bv*10.**9.,A))[0]*eV
            fD2=-quad(IntegrandD,0,(6.*np.pi**2./vs)**(1./3.)*c, args=(h,kb,v,c,T,v0,Bv*10.**9.))[0]*eV

            if j==0: c=np.float64(np.sqrt(Gv0*10.**9./rho0))#
            else: 
                c=np.float64(np.sqrt(Gv*10.**9./rho))
                #c=C*np.sqrt(Bv*10.**9./rho)
                #c=C*np.sqrt(Gv0*Bv/Bv0*10.**9./rho)
                
            y = lambda x: (A**3./(2.*np.pi**2.)*( 1./3.*np.exp(3.*x*(vs)**(1./3.)/(c*A)) - np.exp(2.*x*(vs)**(1./3.)/(c*A)) + np.exp(x*(vs)**(1./3.)/(c*A)) - 1./3. ))-1.
            omega = brentq(y, 10**10., 10**15.)*delta_theta[k]
            theta = omega*h/kb
            
            
            ft3=-quad(Integrand,0,omega, args=(h,kb,c,T,vs,v0,Bv*10.**9.,A))[0]*eV
            fD3=-quad(IntegrandD,0,(6.*np.pi**2./vs)**(1./3.)*c*0.9, args=(h,kb,v,c,T,v0,Bv*10.**9.))[0]*eV

            fs= -(debye_function(T_Deb/T) - 3.*np.log(1.-np.exp(-T_Deb/T)) ) * kb * T + 9./8.*kb*T_Deb



            f.append(1.*ft1+1.*ft2+1.*ft3+9./8.*h*omega*eV+0.0067)

            fd.append(1.*fD1+1.*fD2+1.*fD3+9./8.*(6.*np.pi**2./vs)**(1./3.)*c*h*eV)

            fstd.append(fs*eV)



            F_v.append(f[i])
            FD_v.append(fd[i])
            Fstd_v.append(fstd[i])


            F.append(eos.fbirch(v*10.**(30.),par[0],par[1],par[2],par[3])+f[i])
            FD.append(eos.fbirch(v*10.**(30.),par[0],par[1],par[2],par[3])+fd[i])
            Fstd.append(eos.fbirch(v*10.**(30.),par[0],par[1],par[2],par[3])+fstd[i])
            i+=1



        eF = EOS([v*10.**30. for v in Vexp],F)
        eF.fit_birch()
        pF = eF.par
        vmin=pF[0]*10.**(-30.)
        
        eFd = EOS([v*10.**30. for v in Vexp],FD)
        eFd.fit_birch()
        pFd = eFd.par
        vminD=pFd[0]*10.**(-30.)
        eFstd = EOS([v*10.**30. for v in Vexp],Fstd)
        eFstd.fit_birch()
        pFstd = eFstd.par
        vminstd=pFstd[0]*10.**(-30.)
        print vmin,vminD,vminstd
        Lmin.append(float(vmin*10**30.))
        LminD.append(float(vminD*10**30.))
        Lminstd.append(float(vminstd*10**30.))
        xp=np.linspace(min(V0),max(V0),10)
        #ax2.plot(Vexp,F_v,'b')
        #ax2.plot(Vexp,Fstd_v,'g')
        #plt.plot(V0,f)
        #plt.plot()
        v00=vmin
        rho00=m/vmin
        
        if k==0: f_out.write('{0} {1}\n'.format(T, vmin))
        
        j+=1
        
    #ax1.plot(TEMP,LminD)
    #ax2.plot(TEMP,Lmin,'g')
    trange,al = alpha(list(TEMP),Lmin)
    #ax1.plot(trange,al,'g')

    trange,alD = alpha(list(TEMP),LminD)
    #ax1.plot(trange,alD, 'g--')

    trange,alstd = alpha(list(TEMP),Lminstd)
    #ax1.plot(trange,alstd, 'g-.')
    print al
    ALPHA1.append(al[10])
    ALPHA2.append(al[100])
    ALPHA3.append(al[140])
    ALPHA4.append(al[180])
    ALPHA5.append(al[210])
    ALPHA1e.append(alstd[10])
    ALPHA2e.append(alstd[100])
    ALPHA3e.append(alstd[140])
    ALPHA4e.append(alstd[180])
    ALPHA5e.append(alstd[210])
    ALPHA1a.append(alD[10])
    ALPHA2a.append(alD[100])
    ALPHA3a.append(alD[140])
    ALPHA4a.append(alD[180])
    ALPHA5a.append(alD[210])
    
    

    f_out.close()
    k+=1
    

##############################################################################################

expT = [5,25,50,100,200,293,400,500,600,700,800,900,1000]
expalpha = [0.0006,0.21,0.88,2.6,4.1,4.5,4.5,4.6,4.7,4.8,5.0,5.0,5.2]
#ax1.plot(expT,expalpha,'bo',ms=10., label = r'$W$ $exp.$')


"""
f_Mo = open('Mo_exp.dat')
lines = f_Mo.readlines()
T_Mo = [float(line.split(',')[0]) for line in lines]
alpha_Mo = [float(line.split(',')[1]) for line in lines]
f_Mo.close()


ax1.plot(T_Mo,alpha_Mo,'g', lw=3,label=r'$Mo$ experimental')





expT = [5,25,50,100,200,293,400,500,600,700,800,900,1000]
expalpha = [0.0006,0.21,0.88,2.6,4.1,4.5,4.5,4.6,4.7,4.8,5.0,5.0,5.2]
ax1.plot(expT,expalpha,'bo',ms=10., label = r'$W$ $exp.$')
"""
#ax1.set_ylim(ymax=30.)
#ax1.set_xlabel('T (K)')
#ax2.set_xlabel('T (K)')
#ax1.set_ylabel(r'TEC $(10^{-6} K^{-1})$')
#ax2.set_ylabel(r'$V$ $(\AA^3)$')
con=[0.,3.,6.,9.,12.,18.,25.,50.]
ax1.plot(con,ALPHA1,'b')
ax1.plot(con,ALPHA2,'g')
ax1.plot(con,ALPHA3,'r')
ax1.plot(con,ALPHA4,'c')
ax1.plot(con,ALPHA5,'m')
ax1.plot(con,ALPHA1e,'b-.')
ax1.plot(con,ALPHA2e,'g-.')
ax1.plot(con,ALPHA3e,'r-.')
ax1.plot(con,ALPHA4e,'c-.')
ax1.plot(con,ALPHA5e,'m-.')
ax1.plot(con,ALPHA1a,'b--')
ax1.plot(con,ALPHA2a,'g--')
ax1.plot(con,ALPHA3a,'r--')
ax1.plot(con,ALPHA4a,'c--')
ax1.plot(con,ALPHA5a,'m--')

############
f=open('WRe1_2')
lines=f.readlines()
f.close()

exc=[]
ex400=[]
ex1000=[]
ex1400=[]
ex1800=[]

for l in lines[1::]:
    exc.append(float(l.split(',')[0]))
    ex400.append(float(l.split(',')[1]))
    ex1000.append(float(l.split(',')[2]))
    ex1400.append(float(l.split(',')[3]))
    ex1800.append(float(l.split(',')[4]))
    
ax1.plot(exc,ex400,'go')
ax1.plot(exc,ex1000,'ro')
ax1.plot(exc,ex1400,'co')
ax1.plot(exc,ex1800,'mo')


plt.show()

