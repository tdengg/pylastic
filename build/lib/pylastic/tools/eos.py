import numpy as np
import json
import pickle
import os
import numpy.random
import matplotlib.pyplot as plt
from pylastic.tools import convert_latt_vol
from lxml import etree


class Birch(object):
    def __init__(self, V, E, structure):
        
        
        #%%%%%%%%--- CONSTANTS ---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        _e     = 1.602176565e-19              # elementary charge
        #Bohr   = 5.291772086e-11              # a.u. to meter
        #Ryd2eV = 13.605698066                 # Ryd to eV
        Angstroem = 1.e-10
        self.__cnvrtr = (_e)/(1e9*Angstroem**3)    # Ryd/[a.u.^3] to GPa
        #--------------------------------------------------------------------------------------------------------------------------------
        
        
        self.__E = [e for e in E]
        self.__V = V
        self.__covera=[1.]
        
        self.__mpl = True
        verbose = False
        self.verbose = verbose
        
        self.__structure=structure
        
        ## start parameters:
        self.__b0 = np.float32(60000) # Bulk-Modulus
        self.__db0 = np.float32(4.)   # derivative of Bulk-Modulus with respect to V 
        self.__itmax = 800
        self.__epsilon = 0.00007
        self.__grad = -1
    
    def fit(self):
        
        a = []

        diff = []
        vbad = []
        ebad = []
    
        v0, emin = self.minIn()
        
        arec = np.array([[float(v0), float(self.__b0), float(self.__db0), float(emin)]])
        
        norm_diff1 = np.linalg.norm(arec)
        ressq0, fite0, res0 = (self.fitev(arec, self.__V, self.__E))
        i=0
        while norm_diff1 > self.__epsilon:
            if i<self.__itmax:# aad < 0:
                
                diff0, parnew0 = self.minimization(arec, self.__V, self.__E)
                ressq, fite, res = self.fitev(parnew0, self.__V, self.__E)
                        
                diff1, parnew1 = self.minimization(parnew0, self.__V, self.__E)
                ressq, fite, res = self.fitev(parnew1, self.__V, self.__E)
                arec = parnew1
                        
                norm_diff0 = np.linalg.norm(diff0)
                norm_diff1 = np.linalg.norm(diff1)
                self.__grad = norm_diff1 - norm_diff0
                #print norm_diff1
            else:
                break
            i=i+1
                
        parmin, deltamin = self.minimization2(parnew1, self.__V, self.__E)
        

        #print deltamin
        if deltamin < norm_diff1:
            parnew1 = parmin
            out = '2-Norm of residual vector: ' + str(deltamin)
        else:
            deltamin = norm_diff1
            out = '2-Norm of residual vector: ' + str(deltamin)
        print out
        # convert volume to lattice parameter:
                    
        if self.__structure == 'fcc':
            latt = (4. * parnew1[0,0])**(1./3.)
        elif self.__structure == 'bcc':
            latt = (2. * parnew1[0,0])**(1./3.)
        elif self.__structure == 'hcp':
            covera = input('c over a ratio: ')
            latta = (parnew1[0,0]/(self.__covera*0.866))**(1./3.)
            lattc = latta * self.__covera
        """if verbose == True:
            print('---------------------------------------')
            print('volume:                     ' + str(parnew1[0,0]) + ' Bohr')
            print('---------------------------------------')
            if structure == 'fcc' or structure == 'bcc':
                print('equilibrium lattice parameter:   ' + str(latt*0.529177) + ' Angstroem')
                print('                                 ' + str(latt) + ' Bohr')
            elif structure == 'hcp':
                print('equilibrium lattice parameter a: ' + str(latta*0.529177) + ' Angstroem')
                print('                                 ' + str(latta) + ' Bohr')
                print('                              c: ' + str(lattc*0.529177) + ' Angstroem')
                print('                                 ' + str(lattc) + ' Bohr')
            print('---------------------------------------')
            print('Bulk-Modulus:               ' + str(parnew1[0,1]) + ' au.')
            print('                            ' + str(parnew1[0,1]*2.942104*10**4.) + ' GPa')
            print('---------------------------------------')
            print('derivative of Bulk-Modulus: ' + str(parnew1[0,2]))
            print('---------------------------------------')
            print('minimal energy:             ' + str(parnew1[0,3]) + ' Hartree')
            print('---------------------------------------')"""
        #else:
        #polyfit of volue steps to determine (c/a)min:
        if self.__structure == 'hcp' and len(self.__V) == len(self.__covera):
            coeff = np.polyfit(self.__V, covera, 3)
            poly = np.poly1d(coeff)
            dpoly = np.poly1d.deriv(poly)
            ddpoly = np.poly1d.deriv(dpoly)
            coamin = poly(parnew1[0,0])
            plt.plot(covera,self.__V,coamin,parnew1[0,0])
            #print coamin
        
        convert = convert_latt_vol.Convert(self.__structure)
        
        a0, V0 = convert.volumeToLatt([parnew1[0,0]],self.__covera)
        
        out = (('V0: ' + str(round(parnew1[0,0], 4)))  + ('a0: ' + str(round(a0[0], 4)).rjust(16)+ ('B0: ' + str(round(parnew1[0,1]*self.__cnvrtr, 4))).rjust(16) + ("B0': " + str(round(parnew1[0,2],4))).rjust(16) + ('E0: ' + str(round(parnew1[0,3], 4))).rjust(16)))
        print out
        #if structure == 'fcc':
        #    a0 = (4.*parnew1[0,0])**(1./3.)
        #    out = (('V0: ' + str(round(parnew1[0,0], 4)))  + ('a0: ' + str(round(a0, 4)).rjust(16)+ ('B0: ' + str(round(parnew1[0,1]*2.942104*10**4., 4))).rjust(16) + ("B0': " + str(round(parnew1[0,2],4))).rjust(16) + ('E0: ' + str(round(parnew1[0,3], 4))).rjust(16)))
        #    print out
        #elif structure == 'hcp':
        #    a0 = (2.*parnew1[0,0]/(3.**(1./2.)*float(coamin)))**(1./3.)
        #    out = (('V0: ' + str(round(parnew1[0,0], 4)))  + ('a0: ' + str(a0).rjust(16)+ ('B0: ' + str(round(parnew1[0,1]*2.942104*10**4., 4))).rjust(16) + ("B0': " + str(round(parnew1[0,2],4))).rjust(16) + ('E0: ' + str(round(parnew1[0,3], 4))).rjust(16)))
        #    print out
            
        ##elif structure == 'hex':
        ##    print(('V0: ' + str(round(parnew1[0,0], 4)))  + ('a0: ' + str((2.*parnew1[0,0]/(3.**(1./2.)*float(covera[i])))**(1./3.)).rjust(16)+ ('B0: ' + str(round(parnew1[0,1]*2.942104*10**4., 4))).rjust(16) + ("B0': " + str(round(parnew1[0,2],4))).rjust(16) + ('E0: ' + str(round(parnew1[0,3], 4))).rjust(16)))
        #elif structure == 'bcc':
        #    a0 = (2.*parnew1[0,0])**(1./3.)
        #    out = (('V0: ' + str(round(parnew1[0,0], 4)))  + ('a0: ' + str(round((2*parnew1[0,0])**(1./3.), 4))).rjust(16)+ ('B0: ' + str(round(parnew1[0,1]*2.942104*10**4., 4))).rjust(16) + ("B0': " + str(round(parnew1[0,2],4))).rjust(16) + ('E0: ' + str(round(parnew1[0,3], 4))).rjust(16))
        #    print out
        #elif structure == 'diamond':
        #    a0 = (8.*parnew1[0,0])**(1./3.)
        #    out = (('V0: ' + str(round(parnew1[0,0], 4)))  + ('a0: ' + str(round((8*parnew1[0,0])**(1./3.), 4))).rjust(16)+ ('B0: ' + str(round(parnew1[0,1]*2.942104*10**4., 4))).rjust(16) + ("B0': " + str(round(parnew1[0,2],4))).rjust(16) + ('E0: ' + str(round(parnew1[0,3], 4))).rjust(16))
        #    print out
        #logfile.write(out + '\n')
        #logfile.close()
        #plt.plot(self.__V, fite0)#
        lv = np.linspace(min(self.__V),max(self.__V),100)
        dump, plote, dump = (self.fitev(parnew1, lv, self.__E))
        #print lv, plote
        if self.__mpl:
            plt.cla()
            plt.plot(lv, plote, '')
            plt.plot(self.__V, self.__E, '.')
            plt.xlabel(r'$volume$   $[{Bohr^3}]$')
            plt.ylabel(r'$total$ $energy$   $[{Hartree}]$')
            plt.legend(loc='best')
            self.p = plt
            #plt.savefig(calchome + 'eos.png')
            #plt.show()
            
        #results = etree.Element('plot')
           
        ###############reschild.set(,str(convpar[key]))##############
        self.reschild = etree.Element('graph')
        for i in range(len(lv)):
            point = etree.SubElement(self.reschild, 'point')
            point.set('volume',str(lv[i]))
            point.set('energy',str(plote[i]))
        self.reschild.set('energy_min',str(parnew1[0,3]))
        try:
            self.reschild.set('coa_min',str(coamin))
        except:
            print ''
        self.reschild.set('vol_min',str(parnew1[0,0]))
        self.reschild.set('B0',str(parnew1[0,1]*2.942104*10**4.))
        self.reschild.set('dB0',str(parnew1[0,2]))
        self.reschild2 = etree.Element('graph_exp')
        for i in range(len(self.__V)):
            point2 = etree.SubElement(self.reschild2, 'point')
            point2.set('volume',str(self.__V[i]))
            point2.set('energy',str(self.__E[i]))
        self.reschild3 = etree.Element('graph_exp_bad')
        for i in range(len(vbad)):
            point3 = etree.SubElement(self.reschild3, 'point')
            point3.set('volume',str(vbad[i]))
            point3.set('energy',str(ebad[i]))
        #results.append(reschild2)
        #restree = etree.ElementTree(results)
        #restree.write(calchome + 'eosplot.xml')
        #############################################################
        #plt.show()
                
        self.out0 = parnew1[0,0]
        self.out1 = parnew1[0,1]*2.942104*10**4. #bulk modulus in GPa
        self.out2 = parnew1[0,2]
        self.out3 = parnew1[0,3]
        if len(self.__V) == len(self.__covera):
            self.out4 = coamin
        self.out5 = a0
        self.deltamin = deltamin
        
        self.a = a
        self.v = self.__V
        self.ein = self.__E
        
        if parnew1[0,0] <= min(self.__V)+(max(self.__V)-min(self.__V))*(0.2) or parnew1[0,0] >= max(self.__V)-(max(self.__V)-min(self.__V))*(0.2):
            self.recalculate = True
        else:
            self.recalculate = False
        self.a0 = a0
        
            
    def minIn(self):
        """ find minimum of total energy
        """
        emin = min(self.__E)
        indexv = self.__E.index(emin)
        v0 = np.float32(self.__V[indexv])
        return v0, emin
        
    def fitev(self, par, v, ein):

        fite = []
        deltasq = []
        res = []
        v0 = par[0,0]
        b0 = par[0,1]
        db0 = par[0,2]
        emin = par[0,3]
        i=0
        while i < len(v):
                
            vov = (v0/v[i])**(2./3.)
            fite.append(float(emin + 9. * v0 * b0/16. * ((vov - 1.)**3. * db0 + (vov - 1.)**2. * (6. - 4. * vov))))
            if len(v) == len(ein):
                deltasq.append((fite[i] - ein[i])**2.)
                res.append(fite[i] - ein[i])
            #print (emin - ein[i])**2
            i = i+1
        return deltasq, fite, res
        
    def minimization(self, aold, v, ein):
        i=0
        defit_dV = []
        defit_dB = []
        defit_ddB = []
        defit_demin = []
        jacobian = []
        v0 = aold[0,0]
        b0 = aold[0,1]
        db0 = aold[0,2]
        emin = aold[0,3]
        while i < len(v):
            ## Jacobian: 
            #  derivative of efit with respect to V0
            a = 3. * db0 * ((v0 / v[i]**(2./3.) - v0**(1./3.) )**2. * (1. / v[i]**(2./3.) - 1./3. * v0**(-2./3.)))
            b = 2. * ((v0**(7./6.) / v[i]**(2./3.) - v0**(1./2.)) * (7./6. * v0**(1./6.) / v[i]**(2./3.) - 0.5 * v0 ** (-1./2.)) * (6. - 4. * (v0/v[i])**(2./3.)))
            c = (v0**(7./6.) / v[i]**(2./3.) - v0**(1./2.))**2. * (-8./3. * v0**(-1./3.) / v[i]**(2./3.))
                
            defit_dV.append((-9.)/16.* b0 * (a + b + c))
            #  derivative of efit with respect to B0
            vov = (v0/v[i])**(2./3.)
            defit_dB.append((-9.) * v0 / 16. * ((vov - 1.)**3. * db0 + (vov - 1.)**2. * (6. - 4. * vov)))
            #  derivative of efit with respect to dB0
            defit_ddB.append((-9.) / 16. * b0 * (v0 / v[i]**(2./3.) - v0**(1./3.))**3.)
            #  derivative of efit with respect to emin    
            defit_demin.append(-1)
                
            jacobian.append([defit_dV[i],defit_dB[i],defit_ddB[i],defit_demin[i]])
            ##
            ## residuals:
            ressq, fite, res = self.fitev(aold, v, ein)
                
            i = i+1
        
        
        A = np.matrix(jacobian)
        r = np.array(res)
        B = np.dot(np.transpose(A),A)
        C = np.dot(np.transpose(A),(r))
        delta = np.transpose(np.linalg.solve(B,np.transpose(C)))
        anew = aold + 0.1*delta
        
        return res, anew
        
    def minimization2(self, aold, v, ein):
        delta = []
        anew = np.array([[0.,0.,0.,0.]])
        
        i=0
        while i < 5000:
            anew[0,0] = aold[0,0] * (np.random.uniform(-1,1) * 0.0001 + 1.)
            anew[0,1] = aold[0,1] * (np.random.uniform(-1,1) * 0.0001 + 1.)
            anew[0,2] = aold[0,2] * (np.random.uniform(-1,1) * 0.0001 + 1.)
            anew[0,3] = aold[0,3] * (np.random.uniform(-1,1) * 0.0001 + 1.)
                
            ressq, fite, res = self.fitev(anew, v, ein)
            delta.append(np.linalg.norm(res))
            
            if delta[i]<delta[i-1]:
                deltamin = delta[i]
                amin = anew
                
            else:
                amin = aold
                deltamin = delta[0]
            #amin = anew
            i+=1
        print deltamin  
        return amin, deltamin
    
class Setup(object):
    def __init__(self, Vmin, Vmax, N, cod='vasp', executable='/home/t.dengg/bin/vasp/vasp.5.3/vasp'):
        self.__cod=cod
        if self.__cod=='vasp': from pylastic.io.vasp import POS
        from pylastic.elatoms import Structures, ElAtoms


        ########################## Read in POSCAR file: ######################
        poscar = POS('POSCAR').read_pos()

        ###################### Create Structures instance: ###################
        structures = Structures('vasp')

        ## Generate distorted structures and add them to structures object: ##
        atom = ElAtoms('vasp')
        atom.poscarToAtoms(poscar)
        for scale in np.linspace(Vmin,Vmax,N):

            atom = ElAtoms('vasp')
            atom.poscarToAtoms(poscar)
            atom.scale = scale
            structures.append_structure(atom)
            print atom.scale

        ####################### Write vasp input files: #######################
        structures.write_structures(structures)

        #################### Start local vasp calculation: ####################
        structures.executable = executable
        structures.calc_vasp()

class Analyze(Birch):
    def __init__(self, verbous=False):
        
        from pylastic.postprocess import ECs
        self.__structure='bcc'
        
        
        ec = ECs('vasp')
        ec.set_structures()
        ec.set_gsenergy()
        V=[]
        scale=[]
        gsenergy=[]
        
        for atom in ec.get_atomsByStraintype(None):
            V.append(atom.V)
            scale.append(atom.scale)
            gsenergy.append(atom.gsenergy)
        birch = Birch(V,gsenergy, self.__structure)
        birch.fit()
        plt.plot(V,gsenergy)
        if verbous: plt.show()
        self.__a0 = birch.a0[0]
        
    def set_structure(self, structure):
        self.__structure = structure
        
    def get_structure(self):
        return self.__structure
    
    def get_a0(self):
        return self.__a0
    
    def set_a0(self):
        return self.__a0
    
    a0 = property( fget = get_a0       , fset = set_a0)
    
    