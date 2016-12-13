import os
import re
import sys
import numpy as np

class POS(object):
    """
    EMTO input/output handling.
    
    Note:
    kstr vectors multiplied by SWS_old/SWS_new to give the same volume as before deformation and account for volume change with SWS in kgrn.
    """
    def __init__(self, verbouse=False, OLDKSTR=False):
        self.__OLDKSTR = OLDKSTR
        self.__inputdigits = None
        self.__verbouse=verbouse
        self.__basevect=[]
        self.__lattvect_pos=None
        self.__natoms = 0
        self.__pos={}
        self.__NQ3=1
        self.__N_species=0
        self.__Ry2eV = 13.605698066                 # Ryd to eV
        self.__format = {"kfcd":{"STRNAM":{"read":(10,19) , "write":"{0:9}"},"JOBNAM...":{"read":(10,19) , "write":"{0:9}"}},"shape":{"JOBNAM...":{"read":(10,19) , "write":"{0:9}"}},"kgrn":{"SWS":{"read":(8,14) , "write":"{0:7.4f}"}, "JOBNAM...":{"read":(10,19) , "write":"{0:11}"}}, "kstr":{"A.....":{"read":(10,19) , "write":"{0:10.7f}", "write_nlz":"{0:10.9f}"},"B.....":{"read":(30,39) , "write":"{0:10.7f}", "write_nlz":"{0:10.9f}"},"C.....":{"read":(50,59) , "write":"{0:10.7f}", "write_nlz":"{0:10.9f}"},"Alp":{"read":(10,19) , "write":"{0:10.7f}"},"Bet":{"read":(30,39) , "write":"{0:10.7f}"},"Gam":{"read":(50,59) , "write":"{0:10.7f}"},"BSX":{"read":(10,19) , "write":"{0:10.7f}", "write_nlz":"{0:10.9f}", "write_nnlz":"{0:10.8f}"},"BSY":{"read":(30,39) , "write":"{0:10.7f}", "write_nlz":"{0:10.9f}", "write_nnlz":"{0:10.8f}"},"BSZ":{"read":(50,59) , "write":"{0:10.7f}", "write_nlz":"{0:10.9f}", "write_nnlz":"{0:10.8f}"}, "SWS":{"read":(8,14) , "write":"{0:7.4f}"},"QX":{"read":(10,19) , "write":"{0:10.7f}", "write_nlz":"{0:10.9f}"},"QY":{"read":(30,39) , "write":"{0:10.7f}", "write_nlz":"{0:10.9f}"},"QZ":{"read":(50,59) , "write":"{0:10.7f}", "write_nlz":"{0:10.9f}"},"LAT":{"read":(18,20) , "write":"{0:2.0i}"} , "JOBNAM...":{"read":(10,19) , "write":"{0:9}"}}}
    
    def trans_csystem(self, in_mat, out_mat, in_vec):
        c_vec = np.dot(np.dot(np.linalg.inv(np.transpose(in_mat)),np.transpose(out_mat)), in_vec)
        #c_vec = np.dot(np.dot(np.linalg.inv(in_mat),out_mat), in_vec)
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
        if self.__OLDKSTR==False:
            if parname=='BSX':                                              #New KSTR version
                parname='BS'                                                #New KSTR version
                pos=1
            if parname=='BSY':
                parname='BS'
                pos=2
            if parname=='BSZ':
                parname='BS'
                pos=3                                                       #New KSTR version
            if parname=='QX':                                              #New KSTR version
                parname='QX'                                                #New KSTR version
                pos=1
            if parname=='QY':
                parname='QX'
                pos=2
            if parname=='QZ':
                parname='QX'
                pos=3                                                       #New KSTR version
        for line in ifile:
            if parname in line:
                if parname.startswith('BS') or parname.startswith('QX') and self.__OLDKSTR==False:                            #New KSTR version
                    lstr = line.split()[pos]
                    if self.__inputdigits == None:
                        if lstr[1]=='.':
                            dig=lstr[2:]
                        elif lstr[2]=='.':
                            dig=lstr[3:]
                        else:
                            dig = lstr
                        self.__inputdigits = len(dig)
                    val.append(float(lstr))                #New KSTR version
                else:
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
                if parname in ['QX','BS1','BS2','BS3'] and self.__OLDKSTR==False:
                    if parname == 'QX':
                        for v in val: line = "{0:s} {1:.{4}f} {2:.{5}f} {3:.{6}f}\n".format(parname,v[0],v[1],v[2],str(self.__inputdigits),str(self.__inputdigits),str(self.__inputdigits))
                    else:
                        
                        line = "{0:s} {1:.{4}f} {2:.{5}f} {3:.{6}f}\n".format(parname,val[0],val[1],val[2],str(self.__inputdigits),str(self.__inputdigits),str(self.__inputdigits))
                else:
                    f_read=self.__format[fname][parname]["read"]
                    li = list(line)
                    del li[f_read[0]:f_read[1]]
                    if type(val) is list:
                        if fname=='kstr' and float(val[j]) < 1. and float(val[j]) >= 0.:
                            li[f_read[0]] = self.__format[fname][parname]["write_nlz"].format(float(val[j])).lstrip('0')
                        elif fname=='kstr' and float(val[j]) < 0. and float(val[j]) >= -1.:
                            
                            li[f_read[0]] = '-'+self.__format[fname][parname]["write_nnlz"].format(float(val[j])).lstrip('-0')
                        else:
                            li[f_read[0]] = self.__format[fname][parname]["write"].format(float(val[j]))
                    else: 
                        if fname=='kstr' and float(val) < 1.:
                            li[f_read[0]] = self.__format[fname][parname]["write_nlz"].format(float(val)).lstrip('0')
                        elif fname=='kstr' and float(val) >= 100.: 
                            li[f_read[0]] = self.__format[fname][parname]["write"].format(float(val))[0:-1]
                        else:
                            li[f_read[0]] = self.__format[fname][parname]["write"].format(float(val))
            
                    line = ''.join(li)
                ifile[i] = line
            
                j+=1
            i+=1
        return ifile
    
    def replace_string(self, ifile, fname, parname, val):
        i=0
        j=0
        for line in ifile:
            if parname in line:
        
                f_read=self.__format[fname][parname]["read"]
                li = list(line)
                #del li[f_read[0]:f_read[1]]
            
                
                li[f_read[0]] = self.__format[fname][parname]["write"].format(val)
        
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
        
        #jobnam=self.find(self.__kstr, "kstr", "JOBNAM")
        
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
            self.__pos['latt_mod'] = 'angles'
            BSX, BSY, BSZ = self.calc_basis(latt, A, B, C, alpha, beta, gamma)
        else:
            self.__pos['latt_mod'] = 'vectors'
        
        
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
        #print B_matrix
        if not type(QX)==float:
            for i in range(len(QX)):
                C_vec = np.array([QX[i], QY[i], QZ[i]])
                B_vec.append(self.trans_csystem(B_matrix, C_matrix, C_vec)) #Convert cartesian to direct coordinates
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
            B_vec.append(self.trans_csystem(B_matrix, C_matrix, C_vec))
            self.__pos['vbasis']['b_%s'%(1)] = []
            self.__pos['vbasis']['b_%s'%(1)].append(B_vec[0][0])
            self.__pos['vbasis']['b_%s'%(1)].append(B_vec[0][1])
            self.__pos['vbasis']['b_%s'%(1)].append(B_vec[0][2])
            
            self.__pos['natoms'] = [1]
            
        self.__pos['csystem'] = 'd'
        
        self.__pos['inputdigits'] = self.__inputdigits
        print self.__inputdigits
        
        return self.__pos
        
    
    def set_pos(self,pos):
        self.__pos=pos
    def get_pos(self):
        return self.__pos
    
    def write_kstr(self, pos, fname, posold, pname, STRUCT_name):
        
        self.__inputdigits = posold['inputdigits']
        
        path = fname.rstrip(os.path.basename(fname))
        fname = os.path.basename(fname)
        
        os.system('mkdir %s/%s'%(path,pname))
        path=path+'/%s/'%pname
        os.system('cp kstr.dat %s'%(path+fname))
        kstr=open((path+fname),'rw')
        
        self.__kstr=kstr.readlines()
        
        kstr.close()
        
        #self.__kstr = self.replace_string(self.__kstr, "kstr", "JOBNAM...", STRUCT_name)
        
        
        self.__pos=pos
        BS1 = self.__pos['vlatt_1']
        BS2 = self.__pos['vlatt_2']
        BS3 = self.__pos['vlatt_3']
        #print self.__pos['vlatt_1'],self.__pos['vlatt_2'],self.__pos['vlatt_3']
        B_matrix = np.array([BS1,BS2,BS3])
        C_matrix = np.identity(3)
        
        a=np.sqrt(BS1[0]**2.+BS1[1]**2.+BS1[2]**2.)
        b=np.sqrt(BS2[0]**2.+BS2[1]**2.+BS2[2]**2.)
        c=np.sqrt(BS3[0]**2.+BS3[1]**2.+BS3[2]**2.)
        
        aold=np.sqrt(posold['vlatt_1'][0]**2.+posold['vlatt_1'][1]**2.+posold['vlatt_1'][2]**2.)
        
        #A=a/aold
        #A_scale = self.__pos['vlatt_1'][0]/posold['vlatt_1'][0]
        
        V0latt = abs(np.linalg.det(np.array([posold['vlatt_1'],posold['vlatt_2'],posold['vlatt_3']])))
        cell = np.array([BS1,BS2,BS3])
        self.__V1cell = abs(np.linalg.det(cell))
        A = (V0latt/self.__V1cell)**(1./3.) #Rescale kstr to constant volume
        
        QX=[]
        QY=[]
        QZ=[]
        
        C_vec = []
        for i in range(len(self.__pos['natoms'])):
            B_vec = np.array(self.__pos['vbasis']['b_%s'%(i+1)])
            C_vec.append(self.trans_csystem(C_matrix, B_matrix, B_vec)) #Convert direct to cartesian coordinates
        
        for vec in C_vec:
            QX.append(vec[0])   #everything multiplied by A in previous version
            QY.append(vec[1])   #
            QZ.append(vec[2])   #
        
        if self.__OLDKSTR:
            self.__kstr = self.replace(self.__kstr, "kstr", "QX", QX)                #Old KSTR input
            self.__kstr = self.replace(self.__kstr, "kstr", "QY", QY)                #Old KSTR input
            self.__kstr = self.replace(self.__kstr, "kstr", "QZ", QZ)                #Old KSTR input
        else:
            self.__kstr = self.replace(self.__kstr, "kstr", "QX", C_vec)

        
        
        
        
        alpha = np.arccos(np.dot(BS2,BS3)/(b*c))*180./np.pi
        beta =  np.arccos(np.dot(BS3,BS1)/(a*c))*180./np.pi
        gamma = np.arccos(np.dot(BS1,BS2)/(a*b))*180./np.pi
        
        
        if self.__pos['latt_mod'] == 'angles':
            self.__kstr = self.replace(self.__kstr, "kstr", "A.....", a/a)
            self.__kstr = self.replace(self.__kstr, "kstr", "B.....", b/a)
            self.__kstr = self.replace(self.__kstr, "kstr", "C.....", c/a)
            self.__kstr = self.replace(self.__kstr, "kstr", "Alp", alpha)
            self.__kstr = self.replace(self.__kstr, "kstr", "Bet", beta)
            self.__kstr = self.replace(self.__kstr, "kstr", "Gam", gamma)
        else:
            BSX=[BS1[0],BS2[0],BS3[0]]  #everything multiplied by A in previous version
            BSY=[BS1[1],BS2[1],BS3[1]]  #
            BSZ=[BS1[2],BS2[2],BS3[2]]  #
            if self.__OLDKSTR:
                self.__kstr = self.replace(self.__kstr, "kstr", "BSX", BSX)
                self.__kstr = self.replace(self.__kstr, "kstr", "BSY", BSY)
                self.__kstr = self.replace(self.__kstr, "kstr", "BSZ", BSZ)
            else:
                self.__kstr = self.replace(self.__kstr, "kstr", "BS1", BS1)
                self.__kstr = self.replace(self.__kstr, "kstr", "BS2", BS2)
                self.__kstr = self.replace(self.__kstr, "kstr", "BS3", BS3)
        #self.__kstr = self.replace(self.__kstr,"kstr","BSX", BSX)
        #self.__kstr = self.replace(self.__kstr,"kstr","BSY", BSY)
        #self.__kstr = self.replace(self.__kstr,"kstr","BSZ", BSZ)
        #print self.__kstr
        f=open(path+fname,'w')
        f.writelines(self.__kstr)
        f.close()
        #print a,b,c
        
        self.__natoms = len(QX)
        self.__a = a
        
        
        
        
        V0_Vlatt = (4./3.*self.__natoms*np.pi*posold['scale']**3.)/V0latt
        #print 'V0_latt, V0_sws, V1_cell ',V0latt,(4./3.*self.__natoms*np.pi*posold['scale']**3.), self.__V1cell
        #SWS = self.calc_sws(V0_Vlatt, self.__natoms, self.__V1cell)
        SWS = posold['scale']*(self.__V1cell/V0latt)**(1./3.)
        #print 'SWS ', SWS, posold['scale'], V0latt, np.linalg.det(np.array([BSX,BSY,BSZ]))
        # print a, b, c
        return SWS#a
    
    def read_shape(self, fname):
        if fname: shape = open(fname)
        self.__shape = shape.readlines()
        
        
        return
    
    def write_shape(self, pos, fname, pname, STRUCT_name):
        path = fname.rstrip(os.path.basename(fname))
        fname = os.path.basename(fname)
        
        os.system('mkdir %s/%s'%(path,pname))
        path=path+'/%s/'%pname
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
        
        
        self.__sws = SWS
        self.__pos['scale'] = SWS
        return self.__pos
    
    def calc_volume(self, nat, SWS):
        return 4./3.*nat*np.pi*SWS**3.
    
    def calc_sws(self,V0_Vlatt,nat,V1):
        return (3./(4.*np.pi*nat)*V0_Vlatt*V1)**(1./3.)
    
    def write_kgrn(self, pos, fname, pname, SYSTEM_name):
        
        path = fname.rstrip(os.path.basename(fname))
        fname = os.path.basename(fname)
        
        os.system('mkdir %s/%s'%(path, pname))
        path=path+'/%s/'%pname
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
    
    def write_kfcd(self, pos, fname, pname, SYSTEM_name):
        path = fname.rstrip(os.path.basename(fname))
        fname = os.path.basename(fname)
        
        os.system('mkdir %s/%s'%(path, pname))
        path=path+'/%s/'%pname
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
    def __init__(self, fname = 'etot.dat', funct='GGA'):
        self.__funct = funct
        self.__fname = fname
        self.__gsenergy = 0.
        self.__Ry2eV = 13.605698066                 # Ryd to eV
    
    def set_functional(self,funct):
        if funct in ['LDA','GGA','PBEsol']: self.__funct = funct
        else: print "Sorry, wronc functional%s. Set to 'LDA','GGA' or 'PBEsol' instead"%funct
    
    def get_functional(self):
        return self.__funct
    
    def set_T(self, T):
        self.__T = T
        
    def get_T(self):
        return self.__T

    def set_gsenergy(self):
        """Get groundstate energy from vasprun.xml"""
        if self.__funct == 'LDA': fc = 'LDA'
        elif self.__funct == 'GGA': fc = 'GGA'
        elif self.__funct == 'PBEsol': fc = 'sol'
        #for line in open(self.__fname,'r'): self.__gsenergy = float(line.split()[ind])
        f=open(self.__fname)
        lines = f.readlines()
        f.close()
        for line in lines:
            if 'TOT-%s'%fc in line.split(): self.__gsenergy = float(line.split()[3])#*self.__Ry2eV
        
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
    
    functional = property( fget = get_functional       , fset = set_functional)
    gsenergy = property( fget = get_gsenergy        , fset = set_gsenergy)
    fname = property( fget = get_fname        , fset = set_fname)
    T = property( fget = get_T       , fset = set_T)
    
    