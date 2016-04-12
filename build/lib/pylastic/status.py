import os
import sys
import pickle
import lxml.etree as et

from pylastic.objecthandler import Object
from pylastic.prettyPrint import FileStructure

class Check(FileStructure):
    def __init__(self):
        
        super(FileStructure, self).__init__()
        
        self.__starus = None
        self.__structuresinst = None
        self.__workdir = './'
        self.__cod = 'vasp'
        
        
    
    def set_workdir(self, workdir):
        self.__workdir = workdir
        
    def get_workdir(self):
        return self.__workdir
    
    def set_structuresinst(self, structureinst):
        self.__structuresinst = structureinst
    
    def get_structuresinst(self):
        return self.__structuresinst
    
    def check_calc(self):
        self.__status = {}
        self.__paths = {}
        if not self.__structuresinst:
            if os.path.isfile(self.__workdir+'/structures.pkl'):
                self.__structuresinst = Object().load(self.__workdir+'/structures.pkl')
                atoms = self.__structuresinst.get_structures()
            else:
                raise Exception("Please do setup first!")
                
        else:
            try:
                atoms = self.__structuresinst.get_structures()
            except:
                atoms = self.__structuresinst
        if len(atoms.keys()[0])==2:
            for stype, eta in atoms:
                self.__status[(stype,eta)] = {}
                if self.__cod == 'vasp':
                    try:
                        et.parse(atoms[(stype,eta)].path+'/vasprun.xml')
                        self.__status[(stype,eta)]['status'] = 'finished'
                        atoms[(stype,eta)].status = True
                    except:
                        self.__status[(stype,eta)]['status'] = sys.exc_info()[0]
                        atoms[(stype,eta)].status = False
                
                else: 
                    self.__status[(stype,eta)]['status'] = 'finished'
                    atoms[(stype,eta)].status = True
                self.__status[(stype,eta)]['path'] = atoms[(stype,eta)].path
        
        elif len(atoms.keys()[0])==3:
            for stype, eta, vol in atoms:
                self.__status[(stype,eta,vol)] = {}
                if self.__cod == 'vasp':
                    try:
                        et.parse(atoms[(stype,eta,vol)].path+'/vasprun.xml')
                        self.__status[(stype,eta,vol)]['status'] = 'finished'
                        atoms[(stype,eta,vol)].status = True
                    except:
                        self.__status[(stype,eta,vol)]['status'] = '--------'
                        atoms[(stype,eta,vol)].status = False
                
                else: 
                    self.__status[(stype,eta,vol)]['status'] = 'finished'
                    atoms[(stype,eta,vol)].status = True
                self.__status[(stype,eta,vol)]['path'] = atoms[(stype,eta,vol)].path        
                
            
        self.__statusstring = FileStructure().dicToTree(self.__status)
        f = open('status', 'w')
        f.write(self.__statusstring)
        f.close()
        return self.__status, self.__structuresinst, self.__statusstring
    
    def get_filestructure(self):
        return self.__filestructure
    
    def set_code(self, code):
        if code in ['vasp','exciting','espresso','wien','emto']:
            self.__cod = code
        else:
            print "Unknown code '%s'. Please choose either espresso, exciting, wien, emto or vasp"%code
            
    def get_code(self):
        return self.__cod
    
    workdir = property( fget = get_workdir        , fset = set_workdir)
    structures = property( fget = get_structuresinst       , fset = set_structuresinst)
    code = property( fget = get_code        , fset = set_code    )
    