import os
import sys
import numpy as np
import math


class POS(object):
    
    def __init__(self, case_struct, ordr=2):
        self.__ordr = ordr
        self.__case_struct = case_struct
        pass
    
    def read_in(self):
        
        #%!%!%--- Reading the case.struct file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        INF   = open(self.__case_struct,'r')
        lines = INF.readlines()
        title = lines[0][ 0:80].strip()
        Ltype = lines[1][ 0: 3].strip()
        mode  = lines[2][13:17].strip()
        unit  = lines[2][23:27].strip()
        NAT_L = int(lines[1][27:30])
        
        if (Ltype!='P'   and \
            Ltype!='B'   and \
            Ltype!='F'   and \
            Ltype!='H'   and \
            Ltype!='R'   and \
            Ltype!='CXY' and \
            Ltype!='CXZ'):
            sys.exit('\n.... Oops ERROR: WRONG lattice type, Check the "'+ self.__case_struct +'" file !?!?!?\n')
        
        a1    = float(lines[3][ 0:10])
        a2    = float(lines[3][10:20])
        a3    = float(lines[3][20:30])
        alpha = float(lines[3][30:40])
        beta  = float(lines[3][40:50])
        gamma = float(lines[3][50:60])
        
        if (Ltype == 'H'): 
            gamma  = 120.0
            Ltype  = 'P'
         
        if (Ltype == 'R'):
            a_R    = np.sqrt((a3**2.)/9. + (a1**2.)/3.)
            alpha_R= 2.*math.asin(a1/(2.*a_R))
            a1     = a_R
            a2     = a_R
            a3     = a_R
            alpha  = np.degrees(alpha_R)
            beta   = np.degrees(alpha_R)
            gamma  = np.degrees(alpha_R)
            Ltype  = 'P'
        
        dummy1 = []
        for line in lines:
            if (line.find('MULT=') >= 0): 
                l     = line
                MULT  = int(l[15:17])
                ISPLIT= int(l[34:36])
                dummy1.append([MULT, ISPLIT])
        
        dummy2 = []
        for line in lines:
            if (line.find('NPT=') >= 0):
                l     = line
                Aname = l[0:10]
                NPT   = float(l[15:20])
                R0    = float(l[25:35])
                RMT   = float(l[40:50])
                Z     = float(l[55:60])
                dummy2.append([Aname, NPT, R0, RMT, Z])
        
        dummy3 = []
        nl     = -1
        for line in lines:
            nl += 1
            if (line.find('LOCAL ROT MATRIX:    ') >= 0):
                rotm1 = lines[nl+0][20:41] + lines[nl+0][41:51].strip()
                rotm2 = lines[nl+1][20:41] + lines[nl+1][41:51].strip()
                rotm3 = lines[nl+2][20:41] + lines[nl+2][41:51].strip()
                dummy3.append([rotm1, rotm2, rotm3])
        
        INFO_L = []
        for i in range(NAT_L):
            INFO_L.append(dummy1[i] + dummy2[i] + dummy3[i])
        
        POSN_L = []
        for line in lines:
            if (line.find('X=') >= 0):
                l = line
                X = float(l[12:22])
                Y = float(l[25:35])
                Z = float(l[38:48])
                POSN_L.append([X, Y, Z])
        
        INF.close()
        #--------------------------------------------------------------------------------------------------
        
        
        LC_Dic = {              \
        'CI' :'Cubic I'        ,\
        'CII':'Cubic II'       ,\
        'HI' :'Hexagonal I'    ,\
        'HII':'Hexagonal II'   ,\
        'RI' :'Rhombohedral I' ,\
        'RII':'Rhombohedral II',\
        'TI' :'Tetragonal I'   ,\
        'TII':'Tetragonal II'  ,\
        'O'  :'Orthorhombic'   ,\
        'M'  :'Monoclinic'     ,\
        'N'  :'Triclinic'}
                
        
        path = os.getcwd()
        dir_name   = os.path.basename(path)
        #%!%!%--- Calculating the Space-Group Number and classifying it ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
        os.system('x sgroup')
        sgroup_out = dir_name+'.outputsgroup'
        
        SGf   = open(sgroup_out, 'r')
        SGlins= SGf.readlines()
        SGf.close()
        
        if (SGlins[0].find('warning') >= 0): 
            sys.exit('\n.... Oops WARNING: There is a warning in "'+ sgroup_out +'" file !?!?!?\n')
        
        for i in range(len(SGlins)):
            if (SGlins[i].find('Number and name of space group:') >= 0):
                SGN = int(float(SGlins[i].split()[6]))
                SGN_explanation=SGlins[i].strip()
                break
        
        if (1 <= SGN and SGN <= 2):      # Triclinic
            LC = 'N'
            if (self.__ordr == 2): ECs = 21
            if (self.__ordr == 3): ECs = 56
        
        elif(3 <= SGN and SGN <= 15):    # Monoclinic
            LC = 'M'
            if (self.__ordr == 2): ECs = 13
            if (self.__ordr == 3): ECs = 32
        
        elif(16 <= SGN and SGN <= 74):   # Orthorhombic
            LC = 'O'
            if (self.__ordr == 2): ECs =  9
            if (self.__ordr == 3): ECs = 20
        
        elif(75 <= SGN and SGN <= 88):   # Tetragonal II
            LC = 'TII'
            if (self.__ordr == 2): ECs =  7
            if (self.__ordr == 3): ECs = 16
          
        elif(89 <= SGN and SGN <= 142):  # Tetragonal I
            LC = 'TI'
            if (self.__ordr == 2): ECs =  6
            if (self.__ordr == 3): ECs = 12
        
        elif(143 <= SGN and SGN <= 148): # Rhombohedral II 
            LC = 'RII'
            if (self.__ordr == 2): ECs =  7
            if (self.__ordr == 3): ECs = 20
        
        elif(149 <= SGN and SGN <= 167): # Rhombohedral I
            LC = 'RI'
            if (self.__ordr == 2): ECs =  6
            if (self.__ordr == 3): ECs = 14
        
        elif(168 <= SGN and SGN <= 176): # Hexagonal II
            LC = 'HII'
            if (self.__ordr == 2): ECs =  5
            if (self.__ordr == 3): ECs = 12
        
        elif(177 <= SGN and SGN <= 194): # Hexagonal I
            LC = 'HI'
            if (self.__ordr == 2): ECs =  5
            if (self.__ordr == 3): ECs = 10
        
        elif(195 <= SGN and SGN <= 206): # Cubic II
            LC = 'CII'
            if (self.__ordr == 2): ECs =  3
            if (self.__ordr == 3): ECs =  8
        
        elif(207 <= SGN and SGN <= 230): # Cubic I
            LC = 'CI'
            if (self.__ordr == 2): ECs =  3
            if (self.__ordr == 3): ECs =  6
        else: sys.exit('\n.... Oops ERROR: WRONG Space-Group Number !?!?!?    \n')
        
        if (self.__ordr == 2): order = 'second'
        if (self.__ordr == 3): order = 'third'
        print '\n     '+ SGN_explanation +'\
               \n     '+ LC_Dic[LC] +' structure in the Laue classification.\
               \n     This structure has '+ str(ECs) +' independent '+ order +'-order elastic constants.'
        
        os.system('rm -f :log *outputsgroup* *struct_sgroup')
        #--------------------------------------------------------------------------------------------------

        
        
        #%!%!%--- Preparing the atomic positions of the primitive cell ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        if (Ltype == 'P'  ): NAT_P = NAT_L*1
        if (Ltype == 'B'  ): NAT_P = NAT_L*2
        if (Ltype == 'F'  ): NAT_P = NAT_L*4
        if (Ltype == 'CXY'): NAT_P = NAT_L*2
        if (Ltype == 'CXZ'): NAT_P = NAT_L*2
        
        bgn    = 0
        infl   = []
        INFO_P = []
        POSN_P = []
        for i in range(NAT_L):
            infl = INFO_L[i]
            numa = INFO_L[i][0]
            end  = bgn + numa
            posl = []
            for j in range(bgn, end):
                posl.append(POSN_L[j])
            bgn = end
        
            if (Ltype == 'P'):
                T = [[.0, .0, .0]]
            if (Ltype == 'B'):
                T = [[.0, .0, .0],
                     [.5, .5, .5]]
            if (Ltype == 'F'):
                T = [[.0, .0, .0],
                     [.5, .5, .0],
                     [.5, .0, .5],
                     [.0, .5, .5]]
            if (Ltype == 'CXY'):
                T = [[.0, .0, .0],
                     [.5, .5, .0]]
            if (Ltype == 'CXZ'):
                T = [[.0, .0, .0],
                     [.5, .0, .5]]
        
            for j in range(len(T)):
                INFO_P.append(infl)
                for k in range(len(posl)):
                    POSN_P.append([(T[j][0]+posl[k][0])%1, (T[j][1]+posl[k][1])%1, (T[j][2]+posl[k][2])%1])
        #--------------------------------------------------------------------------------------------------
        
        #%!%!%--- Writing the case_P.struct file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        fP = open(dir_name+'_P.struct', 'w')
        print>>fP, title
        print>>fP, 'P   LATTICE,NONEQUIV.ATOMS:'+'%3d'%(NAT_P)
        print>>fP, 'MODE OF CALC='+mode+' unit='+unit
        print>>fP, '%10.6f'%(a1)   + \
                   '%10.6f'%(a2)   + \
                   '%10.6f'%(a3)   + \
                   '%10.6f'%(alpha)+ \
                   '%10.6f'%(beta) + \
                   '%10.6f'%(gamma)
        
        
        index = '-'
        
        bgn = 0
        for i in range(len(INFO_P)):
            if ( 1 <= i+1 and i+1 <=   9): atom_index ='  '+index+str(i+1)
            if (10 <= i+1 and i+1 <=  99): atom_index = ' '+index+str(i+1)
            if (99 <= i+1 and i+1 <= 999): atom_index =     index+str(i+1)
        
            numa= INFO_P[i][0]
            end = bgn + numa
            for j in range(bgn, end):
                print>>fP, 'ATOM'+atom_index+': X='+'%10.8f'%(POSN_P[j][0]) + \
                                              ' Y='+'%10.8f'%(POSN_P[j][1]) + \
                                              ' Z='+'%10.8f'%(POSN_P[j][2])
                if (j == bgn): print>>fP, '          MULT='  +'%2d'%(INFO_P[i][0]) + \
                                          '          ISPLIT='+'%2d'%(INFO_P[i][1])
            bgn = end
        
            print>>fP, INFO_P[i][2]+' NPT='+'%5d'   %(INFO_P[i][3]) + \
                                    '  R0='+'%10.8f'%(INFO_P[i][4]) + \
                                    ' RMT='+'%10.5f'%(INFO_P[i][5]) + \
                                    '   Z:'+'%7.2f' %(INFO_P[i][6])
            print>>fP, 'LOCAL ROT MATRIX:   '+INFO_P[i][7]
            print>>fP, '                    '+INFO_P[i][8]
            print>>fP, '                    '+INFO_P[i][9]
        print>>fP, '   0      NUMBER OF SYMMETRY OPERATIONS'
        fP.close()
        #--------------------------------------------------------------------------------------------------

        #%!%!%--- Making the M_old matrix ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
        alpha = math.radians(alpha)
        beta  = math.radians(beta)
        gamma = math.radians(gamma)
        
        M_old = np.zeros((3,3))
        
        if (LC == 'CI'  or \
            LC == 'CII' or \
            LC == 'TI'  or \
            LC == 'TII' or \
            LC == 'O'):
        
            M_old[0,0] = a1
            M_old[0,1] = 0.
            M_old[0,2] = 0.
        
            M_old[1,0] = 0.
            M_old[1,1] = a2
            M_old[1,2] = 0.
        
            M_old[2,0] = 0.
            M_old[2,1] = 0.
            M_old[2,2] = a3
        
        
        if (LC == 'HI' or \
            LC == 'HII'):
        
            M_old[0,0] = a1
            M_old[0,1] = 0.
            M_old[0,2] = 0.
        
            M_old[1,0] =-a2/2.
            M_old[1,1] = a2*math.sqrt(3.)/2.
            M_old[1,2] = 0.
        
            M_old[2,0] = 0.
            M_old[2,1] = 0.
            M_old[2,2] = a3
        
        
        if (LC == 'RI' or \
            LC == 'RII'):
            if (SGN == 146 or \
                SGN == 148 or \
                SGN == 155 or \
                SGN == 160 or \
                SGN == 161 or \
                SGN == 166 or \
                SGN == 167 ):
        
                M_old[0,0] = a1*math.sin(alpha/2.)
                M_old[0,1] =-a1*math.sin(alpha/2.)/math.sqrt(3.)
                M_old[0,2] = a1*math.sqrt(1.-(4./3.*(math.sin(alpha/2.))**2.))
        
                M_old[1,0] = 0.
                M_old[1,1] =-M_old[0,1]*2.
                M_old[1,2] = M_old[0,2]
        
                M_old[2,0] =-M_old[0,0]
                M_old[2,1] = M_old[0,1]
                M_old[2,2] = M_old[0,2]
        
            else:
                M_old[0,0] = a1
                M_old[0,1] = 0.
                M_old[0,2] = 0.
        
                M_old[1,0] =-a2/2.
                M_old[1,1] = a2*math.sqrt(3.)/2.
                M_old[1,2] = 0.
        
                M_old[2,0] = 0.
                M_old[2,1] = 0.
                M_old[2,2] = a3
        
        
        if (LC == 'M'):
            M_old[0,0] = a1*math.sin(gamma)
            M_old[0,1] = a1*math.cos(gamma)
            M_old[0,2] = 0.
        
            M_old[1,0] = 0.
            M_old[1,1] = a2
            M_old[1,2] = 0.
        
            M_old[2,0] = 0.
            M_old[2,1] = 0.
            M_old[2,2] = a3
        
        
        if (LC == 'N'):
            M_old[0,0] = a1
            M_old[0,1] = 0.
            M_old[0,2] = 0.
        
            M_old[1,0] = a2*math.cos(gamma)
            M_old[1,1] = a2*math.sin(gamma)
            M_old[1,2] = 0.  
        
            M_old[2,0] = a3*math.cos(beta)
            M_old[2,1] =(a3*(math.cos(alpha)-math.cos(beta)*math.cos(gamma)))/math.sin(gamma)    
            M_old[2,2] = a3*math.sqrt(1.-(math.cos(alpha))**2.-(math.cos(beta))**2.-(math.cos(gamma))**2. \
                       + 2.*math.cos(alpha)*math.cos(beta)*math.cos(gamma))/math.sin(gamma)
        
        D = np.linalg.det(M_old)
        if (Ltype == 'P'  ): V0 = D/1.
        if (Ltype == 'B'  ): V0 = D/2.
        if (Ltype == 'F'  ): V0 = D/4.
        if (Ltype == 'CXY'): V0 = D/2.
        if (Ltype == 'CXZ'): V0 = D/2.
        
        
        A1    = math.sqrt(M_old[0,0]**2. + M_old[0,1]**2. + M_old[0,2]**2.)
        A2    = math.sqrt(M_old[1,0]**2. + M_old[1,1]**2. + M_old[1,2]**2.)
        A3    = math.sqrt(M_old[2,0]**2. + M_old[2,1]**2. + M_old[2,2]**2.)

        ALPHA = math.degrees(math.acos((M_old[1,0]*M_old[2,0] +\
                              M_old[1,1]*M_old[2,1] +\
                              M_old[1,2]*M_old[2,2])/(A2*A3)))
        BETA  = math.degrees(math.acos((M_old[0,0]*M_old[2,0] +\
                              M_old[0,1]*M_old[2,1] +\
                              M_old[0,2]*M_old[2,2])/(A1*A3)))
        GAMMA = math.degrees(math.acos((M_old[0,0]*M_old[1,0] +\
                              M_old[0,1]*M_old[1,1] +\
                              M_old[0,2]*M_old[1,2])/(A1*A2)))

        
        i_dict = {}
        i_dict['sgn'] = SGN
        i_dict['name'] = title
        i_dict['path'] = dir_name
        i_dict['scale'] = 1.
        i_dict['a1'] = A1
        i_dict['a2'] = A2
        i_dict['a3'] = A3
        i_dict['alpha'] = ALPHA
        i_dict['beta'] = BETA
        i_dict['gamma'] = GAMMA
        i_dict['vlatt_1'] = [M_old[0,0],M_old[0,1],M_old[0,2]]
        i_dict['vlatt_2'] = [M_old[1,0],M_old[1,1],M_old[1,2]]
        i_dict['vlatt_3'] = [M_old[2,0],M_old[2,1],M_old[2,2]]
        i_dict['natoms'] = 1
        i_dict["selective"] = None
        i_dict["vbasis"] = {}
        
        
        
        
        return
    
    def write_in(self, i_dict, dir_name):
        
        
        
        A1    = math.sqrt(i_dict['vlatt_1'][0]**2. + i_dict['vlatt_1'][1]**2. + i_dict['vlatt_1'][2]**2.)
        A2    = math.sqrt(i_dict['vlatt_2'][0]**2. + i_dict['vlatt_2'][1]**2. + i_dict['vlatt_2'][2]**2.)
        A3    = math.sqrt(i_dict['vlatt_3'][0]**2. + i_dict['vlatt_3'][1]**2. + i_dict['vlatt_3'][2]**2.)

        ALPHA = math.degrees(math.acos((i_dict['vlatt_2'][0]*i_dict['vlatt_3'][0] +\
                              i_dict['vlatt_2'][1]*i_dict['vlatt_3'][1] +\
                              i_dict['vlatt_2'][2]*i_dict['vlatt_3'][2])/(A2*A3)))
        BETA  = math.degrees(math.acos((i_dict['vlatt_1'][0]*i_dict['vlatt_3'][0] +\
                              i_dict['vlatt_1'][1]*i_dict['vlatt_3'][1] +\
                              i_dict['vlatt_1'][2]*i_dict['vlatt_3'][2])/(A1*A3)))
        GAMMA = math.degrees(math.acos((i_dict['vlatt_1'][0]*i_dict['vlatt_2'][0] +\
                              i_dict['vlatt_1'][1]*i_dict['vlatt_2'][1] +\
                              i_dict['vlatt_1'][2]*i_dict['vlatt_2'][2])/(A1*A2)))
        
        fP   = open(dir_name + '_P.struct', 'r')
        Plins= fP.readlines()
        fP.close()
        fo = open('distorted.struct', 'w')

        Plins.pop(0)
        Plins.insert(0, 'distorted \n')

        Plins.pop(3)
        Lattice_Parameters = '%10.6f'%(A1)   +'%10.6f'%(A2)  +'%10.6f'%(A3)\
                           + '%10.6f'%(ALPHA)+'%10.6f'%(BETA)+'%10.6f'%(GAMMA)
        Plins.insert(3, Lattice_Parameters+'\n')

        for i in range(len(Plins)):
            print>>fo, Plins[i],
        fo.close()
        
        return
    
    def read_sgroup(self):
        return
    
    def write_sgroup(self):
        return

class Energy():
    """Get energies from wien2k output file."""
    def __init__(self, fname = 'distorted_Converged.scf'):
        self.__fname = fname

    def set_gsenergy(self):
        for line in open(self.__fname,'r'):
            if (line.find(':ENE  :')>=0):
                self.__gsenergy = float(line.split()[-1])
        
    def get_gsenergy(self):
        return self.__gsenergy
    
    def set_fname(self,fname):
        self.__fname = fname
        
    def get_fname(self):
        return self.__fname
