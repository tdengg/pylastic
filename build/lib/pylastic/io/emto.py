import os
import re
import sys
import numpy as np

class POS(object):
    def __init__(self, verbouse=False):
        self.__verbouse=verbouse
        self.__basevect=[]
        self.__lattvect_pos=None
        self.__natoms = 0
        self.__pos={}
        self.__NQ3=1
        self.__N_species=0
        
        self.__format = {"kgrn":{"SWS":{"read":(8,14) , "write":"{0:7.4f}"}}, "kstr":{"A.....":{"read":(10,19) , "write":"{0:10.7f}"},"B.....":{"read":(30,39) , "write":"{0:10.7f}"},"C.....":{"read":(50,59) , "write":"{0:10.7f}"},"Alp":{"read":(10,19) , "write":"{0:10.7f}"},"Bet":{"read":(30,39) , "write":"{0:10.7f}"},"Gam":{"read":(50,59) , "write":"{0:10.7f}"},"BSX":{"read":(10,19) , "write":"{0:10.7f}"},"BSY":{"read":(30,39) , "write":"{0:10.7f}"},"BSZ":{"read":(50,59) , "write":"{0:10.7f}"}, "SWS":{"read":(8,14) , "write":"{0:7.4f}"},"QX":{"read":(10,19) , "write":"{0:10.7f}"},"QY":{"read":(30,39) , "write":"{0:10.7f}"},"QZ":{"read":(50,59) , "write":"{0:10.7f}"},"LAT":{"read":(18,20) , "write":"{0:2.0i}"} }}
    
    def trans_csystem(self, in_mat, out_mat, in_vec):
        c_vec = np.dot(np.dot(np.linalg.inv(np.transpose(in_mat)),np.transpose(out_mat)), in_vec)
        return c_vec
        
        
    
    def calc_basis(self, latt, A,B,C, ALF,BET,GAM):
        
             
##
##             Simple cubic
##        
        AlF=ALF*np.pi/180.
        BET=BET*np.pi/180.
        GAM=GAM*np.pi/180.
        
        COA=C/A
        BOA=B/A
        BSX=np.zeros(3)
        BSY=np.zeros(3)
        BSZ=np.zeros(3)
        if latt==1:
            BSX[0]=1.
            BSY[0]=0.
            BSZ[0]=0.
            BSX[1]=0.
            BSY[1]=1.
            BSZ[1]=0.
            BSX[2]=0.
            BSY[2]=0.
            BSZ[2]=1.
#
#           Face centred cubic
#
        elif latt==2:
            BSX[0]=0.5
            BSY[0]=0.5
            BSZ[0]=0.0
            BSX[1]=0.0
            BSY[1]=0.5
            BSZ[1]=0.5
            BSX[2]=0.5
            BSY[2]=0.0
            BSZ[2]=0.5
#
#
#           Body centred cubic
#
        elif latt==3:
            BSX[0]=0.5
            BSY[0]=0.5
            BSZ[0]=-0.5
            BSX[1]=-0.5
            BSY[1]=0.5
            BSZ[1]=0.5
            BSX[2]=0.5
            BSY[2]=-0.5
            BSZ[2]=0.5
#
#           Hexagonal
#
        elif latt==4:
            BSX[0]=1.
            BSY[0]=0.
            BSZ[0]=0.
            BSX[1]=-0.5
            BSY[1]=np.sqrt(3.)/2.
            BSZ[1]=0.
            BSX[2]=0.
            BSY[2]=0.
            BSZ[2]=COA

#
#           Simple tetragonal
#
        elif latt==5:
            BSX[0]=1.
            BSY[0]=0.
            BSZ[0]=0.
            BSX[1]=0.
            BSY[1]=1.
            BSZ[1]=0.
            BSX[2]=0.
            BSY[2]=0.
            BSZ[2]=COA
#
#           Body centred tetragonal
#
        elif latt==6:
            BSX[0]=1.0
            BSY[0]=0.0
            BSZ[0]=0.0
            BSX[1]=0.0
            BSY[1]=1.0
            BSZ[1]=0.0
            BSX[2]=0.5
            BSY[2]=0.5
            BSZ[2]=COA/2.
#
#           Trigonal
#
        elif latt==7:
            BSX[0]=0.
            BSY[0]=1.
            BSZ[0]=COA
            BSX[1]=-np.sqrt(3.)/2.
            BSY[1]=-0.5
            BSZ[1]=COA
            BSX[2]=np.sqrt(3.)/2.
            BSY[2]=-0.5
            BSZ[2]=COA
#
#           Simple orthorombic
#
        elif latt==8:
            BSX[0]=1.
            BSY[0]=0.
            BSZ[0]=0.
            BSX[1]=0.
            BSY[1]=BOA
            BSZ[1]=0.
            BSX[2]=0.
            BSY[2]=0.
            BSZ[2]=COA
#
#           Base centered orthorombic
#
        elif latt==9:
            BSX[0]=1./2.
            BSY[0]=-BOA/2.
            BSZ[0]=0.
            BSX[1]=1./2.
            BSY[1]=BOA/2.
            BSZ[1]=0.
            BSX[2]=0.
            BSY[2]=0.
            BSZ[2]=COA
#
#           Body centred orthorombic
#
        elif latt==10:
            BSX[0]=1./2.
            BSY[0]=-BOA/2.
            BSZ[0]=COA/2.
            BSX[1]=1./2.
            BSY[1]=BOA/2.
            BSZ[1]=-COA/2.
            BSX[2]=-1./2.
            BSY[2]=BOA/2.
            BSZ[2]=COA/2.
#
#           Face centred orthorombic
#
        elif latt==11:
            BSX[0]=1./2.
            BSY[0]=0.
            BSZ[0]=COA/2.
            BSX[1]=1./2.
            BSY[1]=BOA/2.
            BSZ[1]=0.
            BSX[2]=0.
            BSY[2]=BOA/2.
            BSZ[2]=COA/2.
#
#           Simple monoclinic
#
        elif latt==12:
            BSX[0]=1.
            BSY[0]=0.
            BSZ[0]=0.
            BSX[1]=BOA*np.cos(GAM)
            BSY[1]=BOA*np.sin(GAM)
            BSZ[1]=0.
            BSX[2]=0.
            BSY[2]=0.
            BSZ[2]=COA
#
#           Base centred monoclinic
#
        elif latt==13:
            BSX[0]=0.
            BSY[0]=-BOA
            BSZ[0]=0.
            BSX[1]=0.5*np.sin(GAM)
            BSY[1]=-0.5*np.cos(GAM)
            BSZ[1]=-0.5*COA
            BSX[2]=0.5*np.sin(GAM)
            BSY[2]=-0.5*np.cos(GAM)
            BSZ[2]=0.5*COA
#
#           Simple triclinic
#
        elif latt==14:
            BSX[0]=1.
            BSY[0]=0.
            BSZ[0]=0.
            BSX[1]=BOA*np.cos(GAM)
            BSY[1]=BOA*np.sin(GAM)
            BSZ[1]=0.
            BSX[2]=COA*np.cos(BET)
            BSY[2]=COA*(np.cos(ALF)-np.cos(BET)*np.cos(GAM))/np.sin(GAM)
            BSZ[2]=COA*np.sqrt((1.-np.cos(GAM)*np.cos(GAM)-np.cos(ALF)*np.cos(ALF)-np.cos(BET)*np.cos(BET)+2.*np.cos(ALF)*np.cos(BET)*np.cos(GAM)))/np.sin(GAM)
        
        return BSX, BSY, BSZ
        
    def find(self, ifile, fname, parname):
        j=0
        val=[]
        for line in ifile:
            if parname in line:
                f_read=self.__format[fname][parname]["read"]
                val.append(float(line[f_read[0]:f_read[1]]))
            
                j+=1
        if len(val)==1: val=val[0]
        
        return val
    
    def replace(self, ifile, fname, parname, val):
        i=0
        j=0
        for line in ifile:
            if parname in line:
        
                f_read=self.__format[fname][parname]["read"]
                li = list(line)
                del li[f_read[0]:f_read[1]]
                if type(val) is list: li[f_read[0]] = self.__format[fname][parname]["write"].format(float(val[j]))
                else: li[f_read[0]] = self.__format[fname][parname]["write"].format(float(val))
        
                line = ''.join(li)
                ifile[i] = line
            
                j+=1
            i+=1
        return ifile
    
    def read_in(self):
        
        return
    
    def read_kstr(self, fname):
        if fname: 
            kstr = open(fname)
            self.__pos['path']=fname
            
        self.__kstr = kstr.readlines()
        A = self.find(self.__kstr, "kstr", "A.....")
        B = self.find(self.__kstr, "kstr", "B.....")
        C = self.find(self.__kstr, "kstr", "C.....")
        
        alpha = self.find(self.__kstr, "kstr", "Alp")
        beta = self.find(self.__kstr, "kstr", "Bet")
        gamma = self.find(self.__kstr, "kstr", "Gam")
        
        BSX = self.find(self.__kstr, "kstr", "BSX")
        BSY = self.find(self.__kstr, "kstr", "BSY")
        BSZ = self.find(self.__kstr, "kstr", "BSZ")
        
        
        QX=self.find(self.__kstr, "kstr", "QX")
        QY=self.find(self.__kstr, "kstr", "QY")
        QZ=self.find(self.__kstr, "kstr", "QZ")
        
        latt = int(self.find(self.__kstr, "kstr", "LAT"))
        
        if BSX==[]:
            BSX, BSY, BSZ = self.calc_basis(latt, A, B, C, alpha, beta, gamma)
        
        
        BS1 = np.array([BSX[0], BSY[0], BSZ[0]])
        BS2 = np.array([BSX[1], BSY[1], BSZ[1]])
        BS3 = np.array([BSX[2], BSY[2], BSZ[2]])
        
        self.__pos['vlatt_1'] = BS1
        self.__pos['vlatt_2'] = BS2
        self.__pos['vlatt_3'] = BS3
        
        self.__pos['scale'] = 1.
        
        
        B_matrix = np.array([BS1,BS2,BS3])
        C_matrix = np.identity(3)
        
        self.__pos['vbasis'] = {}
        B_vec = []
        if not type(QX)==float:
            for i in range(len(QX)):
                C_vec = np.array([QX[i], QY[i], QZ[i]])
                B_vec.append(self.trans_csystem(C_matrix, B_matrix, C_vec)) #Convert cartesian to direct coordinates
                #print C_vec, B_vec
                #print B_matrix, C_matrix, B_vec
            
            
            for i in range(len(QX)):
                
                self.__pos['vbasis']['b_%s'%(i+1)] = []
                self.__pos['vbasis']['b_%s'%(i+1)].append(B_vec[i][0])
                self.__pos['vbasis']['b_%s'%(i+1)].append(B_vec[i][1])
                self.__pos['vbasis']['b_%s'%(i+1)].append(B_vec[i][2])
            
            self.__pos['natoms'] = map(int,list(np.ones(len(QX))))
        else:
            C_vec = np.array([QX, QY, QZ])
            B_vec.append(self.trans_csystem(C_matrix, B_matrix, C_vec))
            self.__pos['vbasis']['b_%s'%(1)] = []
            self.__pos['vbasis']['b_%s'%(1)].append(B_vec[0][0])
            self.__pos['vbasis']['b_%s'%(1)].append(B_vec[0][1])
            self.__pos['vbasis']['b_%s'%(1)].append(B_vec[0][2])
            
            self.__pos['natoms'] = [1]
            
        self.__pos['csystem'] = 'd'
        return self.__pos
        
    
    def set_pos(self,pos):
        self.__pos=pos
    def get_pos(self):
        return self.__pos
    
    def write_kstr(self, pos, fname):
        path = fname.rstrip(os.path.basename(fname))
        fname = os.path.basename(fname)
        
        os.system('mkdir %s/kstr'%(path))
        path=path+'/kstr/'
        os.system('cp kstr.dat %s'%(path+fname))
        kstr=open((path+fname),'rw')
        
        self.__kstr=kstr.readlines()
        
        kstr.close()
        
        self.__pos=pos
        BS1 = self.__pos['vlatt_1']
        BS2 = self.__pos['vlatt_2']
        BS3 = self.__pos['vlatt_3']
        
        B_matrix = np.array([BS1,BS2,BS3])
        C_matrix = np.identity(3)
        
        a=np.sqrt(BS1[0]**2.+BS1[1]**2.+BS1[2]**2.)
        b=np.sqrt(BS2[0]**2.+BS2[1]**2.+BS2[2]**2.)
        c=np.sqrt(BS3[0]**2.+BS3[1]**2.+BS3[2]**2.)
        
        QX=[]
        QY=[]
        QZ=[]
        
        C_vec = []
        for i in range(len(self.__pos['natoms'])):
            B_vec = np.array(self.__pos['vbasis']['b_%s'%(i+1)])
            C_vec.append(self.trans_csystem(B_matrix, C_matrix, B_vec)) #Convert direct to cartesian coordinates
        
        for vec in C_vec:
            QX.append(vec[0]*a)
            QY.append(vec[1]*a)
            QZ.append(vec[2]*a)
        #for i in range(len(QX)):
        self.__kstr = self.replace(self.__kstr, "kstr", "QX", QX)
        self.__kstr = self.replace(self.__kstr, "kstr", "QY", QY)
        self.__kstr = self.replace(self.__kstr, "kstr", "QZ", QZ)

        
        
        
        
        alpha = np.arccos(np.dot(BS2,BS3)/(b*c))*180./np.pi
        beta =  np.arccos(np.dot(BS3,BS1)/(a*c))*180./np.pi
        gamma = np.arccos(np.dot(BS1,BS2)/(a*b))*180./np.pi
        
        self.__kstr = self.replace(self.__kstr, "kstr", "A.....", 1.)
        self.__kstr = self.replace(self.__kstr, "kstr", "B.....", b/a)
        self.__kstr = self.replace(self.__kstr, "kstr", "C.....", c/a)
        
        self.__kstr = self.replace(self.__kstr, "kstr", "Alp", alpha)
        self.__kstr = self.replace(self.__kstr, "kstr", "Bet", beta)
        self.__kstr = self.replace(self.__kstr, "kstr", "Gam", gamma)
        #self.__kstr = self.replace(self.__kstr,"kstr","BSX", BSX)
        #self.__kstr = self.replace(self.__kstr,"kstr","BSY", BSY)
        #self.__kstr = self.replace(self.__kstr,"kstr","BSZ", BSZ)
        #print self.__kstr
        f=open(path+fname,'w')
        f.writelines(self.__kstr)
        f.close()
        return a
    
    def read_shape(self, fname):
        if fname: shape = open(fname)
        self.__shape = shape.readlines()
        
        
        return
    
    def write_shape(self, pos, fname):
        path = fname.rstrip(os.path.basename(fname))
        fname = os.path.basename(fname)
        
        os.system('mkdir %s/shape'%(path))
        path=path+'/shape/'
        os.system('cp shape.dat %s'%(path+fname))
        shape=open((path+fname),'rw')
        
        self.__shape=shape.readlines()
        shape.close()
        f=open(path+fname,'w')
        f.writelines(self.__shape)
        f.close()
        return
    
    def read_kgrn(self, fname):
        if fname: kgrn = open(fname)
        self.__kgrn = kgrn.readlines()
        
        SWS = self.find(self.__kgrn,"kgrn","SWS")
        self.__pos['scale'] = SWS
        return self.__pos
    
    def write_kgrn(self, pos, fname):
        
        path = fname.rstrip(os.path.basename(fname))
        fname = os.path.basename(fname)
        
        os.system('mkdir %s/kgrn'%(path))
        path=path+'/kgrn/'
        os.system('cp kgrn.dat %s'%(path+fname))
        kgrn=open((path+fname),'rw')
        
        self.__kgrn = kgrn.readlines()
        kgrn.close()
        self.__pos=pos
        SWS = self.__pos['scale']
        
        self.__kgrn = self.replace(self.__kgrn,"kgrn","SWS", SWS)
        
        f=open(path+fname,'w')
        f.writelines(self.__kgrn)
        f.close()
        return
    
    def read_kfcd(self, fname):
        
        if fname: self.__kfcd = open(fname)
        lines = self.__kfcd.readlines()
        self.__split_kfcd =  [re.split(r'(\s+)', l) for l in lines]
        return
    
    def write_kfcd(self, pos, fname):
        path = fname.rstrip(os.path.basename(fname))
        fname = os.path.basename(fname)
        
        os.system('mkdir %s/kfcd'%(path))
        path=path+'/kfcd/'
        os.system('cp kfcd.dat %s'%(path+fname))
        kfcd=open((path+fname),'rw')
        
        self.__kfcd = kfcd.readlines()
        kfcd.close()
        
        f=open(path+fname,'w')
        f.writelines(self.__kfcd)
        f.close()

        return
    
    def write_sgroup(self, pos):
        
        f = open('sgroup.in','w')
        
        f.write('P\n')
        
        a = np.sqrt((pos["vlatt_1"][0])**2. + (pos["vlatt_1"][1])**2. + (pos["vlatt_1"][2])**2.) * pos["scale"]
        b = np.sqrt((pos["vlatt_2"][0])**2. + (pos["vlatt_2"][1])**2. + (pos["vlatt_2"][2])**2.) * pos["scale"]
        c = np.sqrt((pos["vlatt_3"][0])**2. + (pos["vlatt_3"][1])**2. + (pos["vlatt_3"][2])**2.) * pos["scale"]
        alpha = np.arccos(np.dot(pos["vlatt_2"]*pos["scale"],pos["vlatt_3"]*pos["scale"])/(b*c))*180./np.pi
        beta = np.arccos(np.dot(pos["vlatt_3"]*pos["scale"],pos["vlatt_1"]*pos["scale"])/(a*c))*180./np.pi
        gamma = np.arccos(np.dot(pos["vlatt_1"]*pos["scale"],pos["vlatt_2"]*pos["scale"])/(a*b))*180./np.pi
        
        f.write(str(a) +' '+ str(b) +' '+ str(c) +' '+ str(alpha) +' '+ str(beta) +' '+ str(gamma) + '\n')
        natom = 0
        for s in pos['natoms']: natom += s
        f.write(str(natom)+'\n')
        if pos["csystem"] in ['d','D']:
            for i in range(len(pos["natoms"])):
                
                    
                f.write(str(pos["vbasis"]["b_" + str(i+1)][0])  + ' ' + str(pos["vbasis"]["b_" + str(i+1)][1])  + ' ' + str(pos["vbasis"]["b_" + str(i+1)][2]) + '\n')
                f.write('Species_' + str(i+1) + '\n')
        else: print 'Basis vectors in Cartesian coordinates not supported yet!!! \n NOT WRITTEN TO sgroup.in!!!!!'
        f.close()
        
    pos = property( fget = get_pos       , fset = set_pos)

        
class Energy(object):
    """Get energies from EMTO output file."""
    def __init__(self, fname = 'etot.dat', funct='PBEsol'):
        self.__funct = funct
        self.__fname = fname
        
    def set_T(self, T):
        self.__T = T
        
    def get_T(self):
        return self.__T

    def set_gsenergy(self):
        """Get groundstate energy from vasprun.xml"""
        if self.__funct == 'LDA': ind = 1
        elif self.__funct == 'GGA': ind = 2
        elif self.__funct == 'PBEsol': ind = 3
        for line in open(self.__fname,'r'): self.__gsenergy = float(line.split()[ind])
        
    def get_gsenergy(self):
        """Return groundstate energy."""
        return self.__gsenergy
    
    def set_phenergy(self, ph_file):
        """Read phonon free energy from phonopy output (default filename: F_TV) for temperature T."""
        g = open(ph_file)
        lines = g.readlines()
        phenergy = []
        T = []
        for temp in lines:
            phenergy.append(float(temp.split()[1])/96.47244)
            T.append(float(temp.split()[0]))
        g.close()
        
        self.__phenergy = phenergy
        self.__T=T
    
    def get_phenergy(self):
        """Return phonon free energy for temperature T"""
        return self.__phenergy
    
    def set_fname(self,fname):
        self.__fname = fname
        
    def get_fname(self):
        return self.__fname
    
    fname = property( fget = get_fname        , fset = set_fname)
    T = property( fget = get_T       , fset = set_T)
    
    