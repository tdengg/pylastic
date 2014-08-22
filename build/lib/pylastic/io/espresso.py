import os
import sys
import shutil
from numpy import *

class POS(object):
    
    def __init__(self, fname='ElaStic_PW.in'):
        self.__fname = fname
        self.__gsenergy = None
    
    def read_in(self):
        
        fi = open(self.__fname, 'r')
        EIF = fi.read()
        fi.close()
        #--------------------------------------------------------------------------------------------------
        
        #%!%!%--- Check whether the essential NAMELISTs exist in the input file ---%!%!%!%!%!%!%!%!%!%!%!%!
        num1 = EIF.upper().find('&CONTROL')
        if (num1 == -1):
            sys.exit('\n.... Oops ERROR: There is NO "&CONTROL" NAMELIST in the input file !?!?!?'\
                     '\n                 Visit http://www.quantum-espresso.org \n')
        num2 = EIF.upper().find('&SYSTEM')
        if (num2 == -1):
            sys.exit('\n.... Oops ERROR: There is NO "&SYSTEM" NAMELIST in the input file !?!?!?'\
                     '\n                 Visit http://www.quantum-espresso.org \n')
        num3 = EIF.upper().find('&ELECTRONS')
        if (num3 == -1):
            sys.exit('\n.... Oops ERROR: There is NO "&ELECTRONS" NAMELIST in the input file !?!?!?'\
                     '\n                 Visit http://www.quantum-espresso.org \n')
        if (num2 <= num1 or num3 <= num2 ):
            sys.exit('\n.... Oops ERROR: NAMELISTs are NOT sequence in the input file !?!?!?'\
                     '\n                 Visit http://www.quantum-espresso.org \n')
        
        num4 = EIF.upper().find('ATOMIC_SPECIES')
        if (num4 == -1):
            sys.exit('\n.... Oops ERROR: "ATOMIC_SPECIES" card must be specified in the input file !?!?!?'\
                     '\n                 Visit http://www.quantum-espresso.org \n')
        num5 = EIF.upper().find('ATOMIC_POSITIONS')
        if (num5 == -1):
            sys.exit('\n.... Oops ERROR: "ATOMIC_POSITIONS" card must be specified in the input file !?!?!?'\
                     '\n                 Visit http://www.quantum-espresso.org \n')
        #--------------------------------------------------------------------------------------------------
        
        #%!%!%--- Reading the title and &CONTROL NAMELIST ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
        title = EIF[0:EIF.upper().find('&CONTROL')].strip()
        
        CONTROL     = EIF[EIF.upper().find('&CONTROL'):EIF.upper().find('&SYSTEM')].strip()
        CONTROL     = CONTROL.replace(',',' ')
        CONTROL     = CONTROL.replace('=',' = ')
        CONTROLlist = CONTROL.strip().split()
        
        if (CONTROL.lower().find('calculation') != -1):
            for i, key in enumerate(CONTROLlist):
                if (key=='calculation'): CONTROLlist[i+2] = "'relax'"
        else:
            CONTROLlist.insert(1,"'relax'")
            CONTROLlist.insert(1,'=')
            CONTROLlist.insert(1,'calculation')
        
        if (CONTROL.lower().find('tstress') != -1 ):
            for i, key in enumerate(CONTROLlist):
                if (key=='tstress'): CONTROLlist[i+2] = '.true.'
        else:
            CONTROLlist.insert(-1,'tstress')
            CONTROLlist.insert(-1,'=')
            CONTROLlist.insert(-1,'.true.')
        
        if (CONTROL.lower().find('tprnfor') != -1):
            for i, key in enumerate(CONTROLlist):
                if (key=='tprnfor'): CONTROLlist[i+2] = '.true.'
        else:
            CONTROLlist.insert(-1,'tprnfor')
            CONTROLlist.insert(-1,'=')
            CONTROLlist.insert(-1,'.true.')
        #--------------------------------------------------------------------------------------------------
        
        #%!%!%--- Reading the &SYSTEM NAMELIST ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        SYSTEM     = EIF[EIF.upper().find('&SYSTEM'):EIF.upper().find('&ELECTRONS')].strip()
        SYSTEM     = SYSTEM.replace(',',' ')
        SYSTEM     = SYSTEM.replace('=',' = ')
        SYSTEMlist = SYSTEM.strip().split()
        
        if (SYSTEM.lower().find('ibrav') == -1):
            sys.exit('\n.... Oops ERROR: "ibrav" must be specified in "&SYSTEM" NAMELIST !?!?!?'\
                     '\n                 Visit http://www.quantum-espresso.org \n')
        for i, key in enumerate(SYSTEMlist):
            if (key == 'ibrav'):
                ibrav = int(SYSTEMlist[i+2])
        
        if (SYSTEM.lower().find('celldm(1)') == -1):
            sys.exit('\n.... Oops ERROR: "celldm(1)" must be specified in "&SYSTEM" NAMELIST !?!?!?'\
                     '\n                 Visit http://www.quantum-espresso.org \n')
        for i, key in enumerate(SYSTEMlist):
            if (key == 'celldm(1)'):
                celldm1 = float(SYSTEMlist[i+2])
        
        #---Free structures--------------------------------------------------------------------------------
        if (ibrav == 0):
            if (EIF.upper().find('CELL_PARAMETERS') == -1 ):
                sys.exit('\n.... Oops ERROR: "CELL_PARAMETERS" must be specified for the case ibrav = 0'\
                         '\n                 Visit http://www.quantum-espresso.org !?!?!? \n')
            for nl, line in enumerate(EIF.split('\n')):
                if (line.find('CELL_PARAMETERS') != -1): break
            if (line.lower().find('hexagonal') != -1):
                sys.exit('\n.... Oops SORRY: "CELL_PARAMETERS" must be specified in "cubic" coordinate.'\
                         '\n                 The "hexagonal" coordinate has NOT been implemented yet.\n')
        
            fi = open(self.__fname,'r')
            for i in range(nl+1):
                fi.readline()
            M  = []
            nl = 0
        
            CP = zeros((3,3))
            while (nl < 3):
                line = fi.readline()
                if (line == ''): break 
                line = line.strip().split()
                if (len(line) == 3):
                    nl +=1
                    M.append(line)
                elif (len(line) == 0): pass
                else: sys.exit('\n.... Oops ERROR: "CELL_PARAMETERS" are NOT defined correctly !?!?!?'\
                               '\n                 Visit http://www.quantum-espresso.org \n')
            
            fi.close()
            if (nl < 3):
                sys.exit('\n.... Oops ERROR: "CELL_PARAMETERS" are NOT complete in the input file !?!?!?'\
                         '\n                 Visit http://www.quantum-espresso.org \n')
            for i in range(3):
                for j in range(3):
                    CP[i,j] = M[i][j]
        
        #-- Cubic structures ------------------------------------------------------------------------------
        if (ibrav == 1 or ibrav == 2 or ibrav == 3):
            celldm2 = 1. 
            celldm3 = 1. 
            celldm4 = 0.
            celldm5 = 0.
            celldm6 = 0.
        
            CP = zeros((3,3))
            if (ibrav == 1): # cP
                CP[0,0] = 1.
                CP[0,1] = 0.
                CP[0,2] = 0.
                CP[1,0] = 0.
                CP[1,1] = 1.
                CP[1,2] = 0.
                CP[2,0] = 0.
                CP[2,1] = 0.
                CP[2,2] = 1.
            if (ibrav == 2): # cF
                CP[0,0] =-.5
                CP[0,1] = 0.
                CP[0,2] = .5
                CP[1,0] = 0.
                CP[1,1] = .5
                CP[1,2] = .5
                CP[2,0] =-.5
                CP[2,1] = .5
                CP[2,2] = 0.
            if (ibrav == 3): # cI
                CP[0,0] = .5
                CP[0,1] = .5
                CP[0,2] = .5
                CP[1,0] =-.5
                CP[1,1] = .5
                CP[1,2] = .5
                CP[2,0] =-.5
                CP[2,1] =-.5
                CP[2,2] = .5
        #---Hexagonal structures---------------------------------------------------------------------------
        if (ibrav == 4):
            if (SYSTEM.lower().find('celldm(3)') == -1):
                sys.exit('\n.... Oops ERROR: "celldm(3)" must be specified for Hexagonal system !?!?!?'\
                         '\n                 Visit http://www.quantum-espresso.org \n')
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(3)'):
                    celldm3 = float(SYSTEMlist[i+2])
            celldm2 = 1. 
            celldm4 = 0.
            celldm5 = 0.
            celldm6 =-.5
        
            CP      = zeros((3,3))
            CP[0,0] = 1.
            CP[0,1] = 0.
            CP[0,2] = 0.
            CP[1,0] =-.5
            CP[1,1] = sqrt(3.)/2.
            CP[1,2] = 0.
            CP[2,0] = 0.
            CP[2,1] = 0.
            CP[2,2] = celldm3
        
        #---Trigonal or Rhombohedragonal structures--------------------------------------------------------
        if (ibrav == 5):
            if (SYSTEM.lower().find('celldm(4)') == -1):
                sys.exit('\n.... Oops ERROR: "celldm(4)" must be specified for Trigonal system !?!?!?'\
                         '\n                 Visit http://www.quantum-espresso.org \n')
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(4)'):
                    celldm4 = float(SYSTEMlist[i+2])
            celldm2 = 1. 
            celldm3 = 1.
            if (celldm4 < -1 or celldm4 > 1.):
                sys.exit('\n.... Oops ERROR: "celldm(4)" is WRONG !?!?!?    \n')
            celldm5 = celldm4 
            celldm6 = celldm4
        
            CP      = zeros((3,3))
            CP[0,0] = sqrt((1.-celldm4)/2.)
            CP[0,1] =-CP[0,0]/sqrt(3.)
            CP[0,2] = sqrt(1.-(2./3.*(1.-celldm4)))
            CP[1,0] = 0.
            CP[1,1] = 2.*CP[0,0]/sqrt(3.)
            CP[1,2] = CP[0,2]
            CP[2,0] =-CP[0,0]
            CP[2,1] = CP[0,1]
            CP[2,2] = CP[0,2]
        
        #---Tetragonal structures--------------------------------------------------------------------------
        if (ibrav == 6 or ibrav == 7):
            if (SYSTEM.lower().find('celldm(3)') == -1):
                sys.exit('\n.... Oops ERROR: "celldm(3)" must be specified for Tetragonal system !?!?!?'\
                         '\n                 Visit http://www.quantum-espresso.org \n')
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(3)'):
                    celldm3 = float(SYSTEMlist[i+2])
            celldm2 = 1.
            celldm4 = .0 
            celldm5 = .0
            celldm6 = .0
        
            CP = zeros((3,3))
            if (ibrav == 6):
                CP[0,0] = 1.
                CP[0,1] = 0.
                CP[0,2] = 0.
                CP[1,0] = 0.
                CP[1,1] = 1.
                CP[1,2] = 0.
                CP[2,0] = 0.
                CP[2,1] = 0.
                CP[2,2] = celldm3
            if (ibrav == 7):
                CP[0,0] = 0.5
                CP[0,1] =-0.5
                CP[0,2] = 0.5*celldm3
                CP[1,0] = 0.5
                CP[1,1] = 0.5
                CP[1,2] = 0.5*celldm3
                CP[2,0] =-0.5
                CP[2,1] =-0.5
                CP[2,2] = 0.5*celldm3
        
        #---Orthorhombic structures------------------------------------------------------------------------
        celldmlist = ['celldm(2)','celldm(3)']
        if (ibrav == 8 or ibrav == 9 or ibrav == 10 or ibrav == 11):
            for celldmi in celldmlist:
                if (SYSTEM.lower().find(celldmi) == -1):
                    sys.exit('\n.... Oops ERROR: "'+ celldmi +'" must be specified for Orthorhombic system !?!?!?'\
                             '\n                 Visit http://www.quantum-espresso.org \n')
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(2)'):
                    celldm2 = float(SYSTEMlist[i+2])
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(3)'):
                    celldm3 = float(SYSTEMlist[i+2])
            celldm4 = .0
            celldm5 = .0
            celldm6 = .0
        
            CP = zeros((3,3))
            if (ibrav == 8):  # Orthorhombic
                CP[0,0] = 1.
                CP[0,1] = 0.
                CP[0,2] = 0.
                CP[1,0] = 0.
                CP[1,1] = celldm2
                CP[1,2] = 0.
                CP[2,0] = 0.
                CP[2,1] = 0.
                CP[2,2] = celldm3
            if (ibrav == 9):  # oC
                CP[0,0] = 0.5
                CP[0,1] = 0.5*celldm2
                CP[0,2] = 0.
                CP[1,0] =-0.5
                CP[1,1] = 0.5*celldm2
                CP[1,2] = 0.
                CP[2,0] = 0.
                CP[2,1] = 0.
                CP[2,2] = celldm3
            if (ibrav == 10):  # oF
                CP[0,0] = 0.5
                CP[0,1] = 0.
                CP[0,2] = 0.5*celldm3
                CP[1,0] = 0.5
                CP[1,1] = 0.5*celldm2
                CP[1,2] = 0.
                CP[2,0] = 0.
                CP[2,1] = 0.5*celldm2
                CP[2,2] = 0.5*celldm3
            if (ibrav == 11):  # oI
                CP[0,0] = 0.5
                CP[0,1] = 0.5*celldm2
                CP[0,2] = 0.5*celldm3
                CP[1,0] =-0.5
                CP[1,1] = 0.5*celldm2
                CP[1,2] = 0.5*celldm3
                CP[2,0] =-0.5
                CP[2,1] =-0.5*celldm2
                CP[2,2] = 0.5*celldm3
        
        #---Monoclinic structures--------------------------------------------------------------------------
        celldmlist = ['celldm(2)','celldm(3)','celldm(4)']
        if (ibrav == 12 or ibrav == 13):
            for celldmi in celldmlist:
                if (SYSTEM.lower().find(celldmi) == -1):
                    sys.exit('\n.... Oops ERROR: "'+celldmi+'" must be specified for Monoclinic system !?!?!?'\
                             '\n                 Visit http://www.quantum-espresso.org \n')
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(2)'):
                    celldm2 = float(SYSTEMlist[i+2])
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(3)'):
                    celldm3 = float(SYSTEMlist[i+2])
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(4)'):
                    celldm4 = float(SYSTEMlist[i+2])
            if (celldm4 < -1 or celldm4 > 1.):
                sys.exit('\n.... Oops ERROR: "celldm(4)" is WRONG !?!?!?    \n')
            celldm5 = .0
            celldm6 = .0
        
            CP = zeros((3,3))
            if (ibrav == 12):  # mP
                CP[0,0] = 1.
                CP[0,1] = 0.
                CP[0,2] = 0.
                CP[1,0] = celldm2*celldm4
                CP[1,1] = celldm2*sqrt(1.-celldm4*celldm4)
                CP[1,2] = 0.
                CP[2,0] = 0.
                CP[2,1] = 0.
                CP[2,2] = celldm3
            if (ibrav == 13):  # mS
                CP[0,0] = 0.5
                CP[0,1] = 0.
                CP[0,2] =-0.5*celldm3
                CP[1,0] = celldm2*celldm4
                CP[1,1] = celldm2*sqrt(1.-celldm4*celldm4)
                CP[1,2] = 0.
                CP[2,0] = 0.5
                CP[2,1] = 0.
                CP[2,2] = 0.5*celldm3
        
        #---Triclinic structures---------------------------------------------------------------------------
        celldmlist = ['celldm(2)','celldm(3)','celldm(4)','celldm(5)','celldm(6)']
        if (ibrav == 14):
            for celldmi in celldmlist:
                if (SYSTEM.lower().find(celldmi) == -1):
                    sys.exit('\n.... Oops ERROR: "'+celldmi+'" must be specified for Triclinic system !?!?!?'\
                             '\n                 Visit http://www.quantum-espresso.org \n')
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(2)'):
                    celldm2 = float(SYSTEMlist[i+2])
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(3)'):
                    celldm3 = float(SYSTEMlist[i+2])
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(4)'):
                    celldm4 = float(SYSTEMlist[i+2])
            if (celldm4 < -1 or celldm4 > 1.):
                sys.exit('\n.... Oops ERROR: "celldm(4)" is WRONG !?!?!?    \n')
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(5)'):
                    celldm5 = float(SYSTEMlist[i+2])
            if (celldm5 < -1 or celldm5 > 1.):
                sys.exit('\n.... Oops ERROR: "celldm(5)" is WRONG !?!?!?    \n')
            for i, key in enumerate(SYSTEMlist):
                if (key == 'celldm(6)'):
                    celldm6 = float(SYSTEMlist[i+2])
            if (celldm6 < -1 or celldm6 > 1.):
                sys.exit('\n.... Oops ERROR: "celldm(6)" is WRONG !?!?!?    \n')
            
            CP      = zeros((3,3))
            CP[0,0] = 1.
            CP[0,1] = 0.
            CP[0,2] = 0.
            CP[1,0] = celldm2*celldm6
            CP[1,1] = celldm2*sqrt(1.-celldm6**2.)
            CP[1,2] = 0.
            CP[2,0] = celldm3*celldm5
            CP[2,1] =(celldm3*(celldm4-celldm5*celldm6))/(sqrt(1.-celldm6**2.))
            CP[2,2] = celldm3*sqrt(1.-celldm4**2.-celldm5**2.-celldm6**2.+\
                      2.*celldm4*celldm5*celldm6)/(sqrt(1.-celldm6**2.))
        
        #--------------------------------------------------------------------------------------------------
        for i, key in enumerate(SYSTEMlist):
            if (key=='ibrav'):
                SYSTEMlist[i+2] = '0'
        
        celldmlist = ['celldm(2)','celldm(3)','celldm(4)','celldm(5)','celldm(6)']
        for celldmi in celldmlist:
            if (SYSTEM.lower().find(celldmi) != -1):
                for i, key in enumerate(SYSTEMlist):
                    if (key==celldmi):
                        SYSTEMlist.pop(i)
                        SYSTEMlist.pop(i)
                        SYSTEMlist.pop(i)
        
        if (SYSTEM.lower().find('nat') == -1):
            sys.exit('\n.... Oops ERROR: "nat" must be specified in "&SYSTEM" NAMELIST !?!?!?'\
                     '\n                 Visit http://www.quantum-espresso.org \n')
        for i, key in enumerate(SYSTEMlist):
            if (key == 'nat'):
                nat = int(SYSTEMlist[i+2])
        
        if (SYSTEM.lower().find('ntyp') == -1):
            sys.exit('\n.... Oops ERROR: "ntyp" must be specified in the "&SYSTEM" NAMELIST !?!?!?'\
                     '\n                 Visit http://www.quantum-espresso.org \n')
        for i, key in enumerate(SYSTEMlist):
            if (key == 'ntyp'):
                ntyp = int(SYSTEMlist[i+2])
        #--------------------------------------------------------------------------------------------------
        
        #%!%!%--- Reading the &ELECTRONS NAMELIST ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
        TMP = EIF[EIF.upper().find('&ELECTRONS')+1:-1]
        TMP = TMP.replace('ATOMIC_SPECIES','&ATOMIC_SPECIES')
        ELECTRONS = TMP[TMP.upper().find('ELECTRONS'):TMP.find('&')].strip()
        ELECTRONS = ELECTRONS.replace(',',' ')
        ELECTRONS = ELECTRONS.replace('=',' = ')
        ELECTRONSlist = ELECTRONS.strip().split()
        #--------------------------------------------------------------------------------------------------
        
        #%!%!%--- Reading the &IONS NAMELIST ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        if (EIF.upper().find('&IONS') != -1):
            TMP = EIF[EIF.upper().find('&IONS')+1:-1]
            TMP = TMP.replace('ATOMIC_SPECIES','&ATOMIC_SPECIES')
            IONS= TMP[TMP.upper().find('IONS'):TMP.find('&')].strip()
        
            IONS     = IONS.replace(',',' ')
            IONS     = IONS.replace('=',' = ')
            IONSlist = IONS.strip().split()
        
        else:    
            IONSlist = ['&IONS','ion_dynamics','=',"'bfgs'",'/']
        #--------------------------------------------------------------------------------------------------
        
        #%!%!%--- Reading the ATOMIC_SPECIES, ATOMIC_POSITIONS, and K_POINT INPUT_CARDS ---%!%!%!%!%!%!%!%!
        ATOMIC_SPECIES = EIF[EIF.upper().find('ATOMIC_SPECIES'):EIF.upper().find('ATOMIC_POSITIONS')-1].strip()
        
        for nl, line in enumerate(EIF.split('\n')):
            if (line.find('ATOMIC_POSITIONS') != -1): break
        if (line.lower().find('crystal') == -1):
            sys.exit('\n.... Oops SORRY: "ATOMIC_POSITIONS" must be specified in "crystal" coordinate'\
                     '\n                 alat, bohr, and angstrom have NOT been implemented yet.\n')
        
        fi = open(self.__fname, 'r')
        for i in range(nl+1):
            fi.readline()
        ATPO = []
        nl   = 0
        while (nl < nat):
            line = fi.readline()
            if (line == ''): break 
            line = line.strip().split()
            if (len(line) == 4):
                nl +=1
                ATPO.append(line)
            elif (len(line) == 0): pass
            else:
                sys.exit('\n.... Oops ERROR: "ATOMIC_POSITIONS" are NOT defined correctly !?!?!?'\
                         '\n                 Visit http://www.quantum-espresso.org \n')
        fi.close()
        if (nl < nat):
            sys.exit('\n.... Oops ERROR: "ATOMIC_POSITIONS" are NOT complete in the input file !?!?!?\
                      \n                 Visit http://www.quantum-espresso.org \n')
        OATPO = []
        while (len(ATPO) > 0):
            natpo = len(ATPO)
            ATYP  = ATPO[natpo-1][0]
            for i in range(natpo-1,-1,-1):
                if (ATYP == ATPO[i][0]):
                    OATPO.append(ATPO[i])
                    ATPO.pop(i)
        
        K_POINTS = EIF[EIF.upper().find('K_POINTS'):EIF.upper().find('CELL_PARAMETERS')].strip()
        #--------------------------------------------------------------------------------------------------
        
        #%!%!%--- Converting the QUANTUM-ESPRESSO input file to exiting input file ---%!%!%!%!%!%!%!%!%!%!%
        #fmt  = '%15.10f'
        #Xi = open('input.xml','w')
        #print >>Xi,'<input>'
        #print >>Xi,'   <structure >'
        #print >>Xi,'      <crystal scale="'+str(celldm1)+'">'
        #for i in range(3):
        #    print >>Xi,'         <basevect>',fmt%CP[i,0],fmt%CP[i,1],fmt%CP[i,2],'</basevect>'
        #print >>Xi,'      </crystal>'
        
        #atom = str(OATPO[nat-1][0])
        #print >>Xi,'      <species speciesfile = "'+atom+'">'
        #for i in range(nat-1,-1,-1):
        #    if (atom == str(OATPO[i][0])):
        #        print >>Xi,'         <atom coord=" ',fmt%(float(OATPO[i][1]))\
        #                                            ,fmt%(float(OATPO[i][2]))\
        #                                            ,fmt%(float(OATPO[i][3])),' "/>'
        #    else:
        #        print >>Xi,'      </species>'
        #        atom = str(OATPO[i][0])
        #        print >>Xi,'      <species speciesfile = "'+atom+'">'
        #        print >>Xi,'         <atom coord=" ',fmt%(float(OATPO[i][1]))\
        #                                            ,fmt%(float(OATPO[i][2]))\
        #                                            ,fmt%(float(OATPO[i][3])),' "/>'
        #print >>Xi,'      </species>'
        #print >>Xi,'   </structure>'
        #print >>Xi,'</input>'
        #Xi.close()
        #os.system("xsltproc $ElaSticROOT/exciting2sgroup.xsl input.xml > sgroup.in;sgroup sgroup.in > sgroup.out")
        #--------------------------------------------------------------------------------------------------
        
        #%!%!%--- Writing the "sgroup.in" file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
        a1 = celldm1 * sqrt(CP[0,0]**2 + CP[0,1]**2 + CP[0,2]**2)
        a2 = celldm1 * sqrt(CP[1,0]**2 + CP[1,1]**2 + CP[1,2]**2)
        a3 = celldm1 * sqrt(CP[2,0]**2 + CP[2,1]**2 + CP[2,2]**2)
        alpha = degrees(math.acos((CP[1,0]*CP[2,0]+CP[1,1]*CP[2,1]+CP[1,2]*CP[2,2])*celldm1**2/(a2*a3)))
        beta  = degrees(math.acos((CP[0,0]*CP[2,0]+CP[0,1]*CP[2,1]+CP[0,2]*CP[2,2])*celldm1**2/(a1*a3)))
        gamma = degrees(math.acos((CP[0,0]*CP[1,0]+CP[0,1]*CP[1,1]+CP[0,2]*CP[1,2])*celldm1**2/(a1*a2)))
        
        si = open('sgroup.in','w')
        print >>si,'P'
        print >>si,a1, a2, a3, alpha, beta, gamma
        print >>si
        print >>si, nat
        for i in range(nat-1,-1,-1):
            print >>si, '%15.10f'%(float(OATPO[i][1])), \
                        '%15.10f'%(float(OATPO[i][2])), \
                        '%15.10f'%(float(OATPO[i][3]))
            print >>si, str(OATPO[i][0])
        si.close()
        
        
        #%!%!%--- Calculating the Space-Group Number and classifying it ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
        os.system('sgroup sgroup.in 1>sgroup.out 2>sgroup.err; rm -f sgroup.in')
        
        if (os.path.getsize('sgroup.err') != 0):
            fer  = open('sgroup.err', 'r')
            lines= fer.readlines()
            print '\n.... Oops '+ lines[0]
            for i in range(1, len(lines)):
                print '                 '+ lines[i]
            fer.close()
            sys.exit()
        else: os.system('rm -f sgroup.err')
        
        SGf   = open('sgroup.out', 'r')
        SGlins= SGf.readlines()
        SGf.close()
        
        for i in range(len(SGlins)):
            if (SGlins[i].find('Number and name of space group:') >= 0):
                SGN = int(float(SGlins[i].split()[6]))
                SGN_explanation = SGlins[i].strip()
                break
          
        i_dict = {}
        i_dict['sgn'] = SGN
        i_dict['name'] = title
        i_dict['path'] = os.getcwd()
        i_dict['scale'] = 1
        i_dict['vlatt_1'] = [CP[0,0],CP[0,1],CP[0,2]]
        i_dict['vlatt_2'] = [CP[1,0],CP[1,1],CP[1,2]]
        i_dict['vlatt_3'] = [CP[2,0],CP[2,1],CP[2,2]]
        i_dict['natoms'] = nat
        i_dict["selective"] = None
        i_dict['parameters'] = {}
        i_dict['parameters']['ELECTRONSlist'] = ELECTRONSlist
        i_dict['parameters']['IONSlist'] = IONSlist
        i_dict['parameters']['SYSTEMlist'] = SYSTEMlist
        i_dict['parameters']['CONTROLlist'] = CONTROLlist
        i_dict['parameters']['ATOMIC_SPECIES'] = ATOMIC_SPECIES
        i_dict['parameters']['K_POINTS'] = K_POINTS
        i_dict['parameters']['OATPO'] = OATPO
        i_dict['parameters']['nat'] = nat
        
        i_dict["vbasis"] = {}
        
        
        
        return i_dict
    
    def write_in(self, i_dict, path):
        
        ELECTRONSlist = i_dict['parameters']['ELECTRONSlist']
        IONSlist = i_dict['parameters']['IONSlist']
        SYSTEMlist = i_dict['parameters']['SYSTEMlist']
        CONTROLlist = i_dict['parameters']['CONTROLlist']
        ATOMIC_SPECIES = i_dict['parameters']['ATOMIC_SPECIES']
        K_POINTS = i_dict['parameters']['K_POINTS']
        nat = i_dict['parameters']['nat']
        OATPO = i_dict['parameters']['OATPO']
        print path
        shutil.copy(os.getcwd()+'/espresso.in',path)
        
        pwi = open(path, 'w')
        
        if (len(i_dict['name']) != 0):
            print >>pwi, i_dict['name']
        
        print >>pwi,'&CONTROL'
        for i in range(1,len(CONTROLlist)-1,3):
            if (len(CONTROLlist[i]) < 17):
                print >>pwi, '   '+'%-16s'%(CONTROLlist[i])+' = '+CONTROLlist[i+2]
            else:
                print >>pwi, '   '+CONTROLlist[i]+' = '+CONTROLlist[i+2]
        print >>pwi,' /'
        
        print >>pwi,'&SYSTEM'
        for i in range(1,len(SYSTEMlist)-1,3):
            if (len(SYSTEMlist[i]) < 17):
                print >>pwi, '   '+'%-16s'%(SYSTEMlist[i])+' = '+SYSTEMlist[i+2]
            else:
                print >>pwi, '   '+SYSTEMlist[i]+' = '+SYSTEMlist[i+2]
        print >>pwi,' /'
        
        print >>pwi,'&ELECTRONS'
        for i in range(1,len(ELECTRONSlist)-1,3):
            if (len(SYSTEMlist[i]) < 17):
                print >>pwi, '   '+'%-16s'%(ELECTRONSlist[i])+' = '+ELECTRONSlist[i+2]
            else:
                print >>pwi, '   '+ELECTRONSlist[i]+' = '+ELECTRONSlist[i+2]
        print >>pwi,' /'
        
        print >>pwi,'&IONS'
        for i in range(1,len(IONSlist)-1,3):
            if (len(IONSlist[i]) < 17):
                print >>pwi, '   '+'%-16s'%(IONSlist[i])+' = '+IONSlist[i+2]
            else:
                print >>pwi, '   '+IONSlist[i]+' = '+IONSlist[i+2]
        print >>pwi,' /'
        
        print >>pwi, ATOMIC_SPECIES
        
        print >>pwi, 'ATOMIC_POSITIONS (crystal)'
        for i in range(nat-1,-1,-1):
            print >>pwi,'   '+'%-3s'%(OATPO[i][0]), '%15.10f'%(float(OATPO[i][1])), '%15.10f'%(float(OATPO[i][2])), \
                                                                                    '%15.10f'%(float(OATPO[i][3]))
        print >>pwi, K_POINTS
        
        print >>pwi, 'CELL_PARAMETERS (cubic)'
        for i in range(3):
            print >>pwi, '%15.10f'%i_dict['vlatt_%s'%str(i+1)][0], '%15.10f'%i_dict['vlatt_%s'%str(i+1)][1], '%15.10f'%i_dict['vlatt_%s'%str(i+1)][2]
        pwi.close()
        
        #fi  = open(os.getcwd()+'/'+path,'r')
        #TMPi= fi.read()
        #TMP = TMPi[TMPi.find('&CONTROL'):TMPi.find('CELL_PARAMETERS')].strip()
        #TMP = TMP + '\nCELL_PARAMETERS (cubic)'
        #fi.close()
        
        return
    def set_gsenergy(self):
        for line in open(self.__fname,'r'):
            if (line.find('!    total energy')>=0):
                self.__gsenergy = float(line.split()[-2])
        
    def get_gsenergy(self):
        return self.__gsenergy
    
    def read_sgroup(self):
        return
    
    def write_sgroup(self):
        return
    
    