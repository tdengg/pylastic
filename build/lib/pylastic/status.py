import os
import sys
import pickle
import lxml.etree as et
from pylastic.objecthandler import Object

class Check(object):
    def __init__(self):
        self.__starus = None
        self.__structuresinst = None
        self.__workdir = './'
        
    def set_structureinst(self, structureinst):
        self.__structuresinst = structureinst
    
    def get_structureinst(self):
        return self.__structuresinst
    
    def check_calc(self):
        self.__status = {}
        self.__paths = {}
        if not self.__structuresinst:
            if os.path.isfile(self.__workdir+'/structures.pkl'):
            #    with open(self.__workdir+'/structures.pkl', 'rb') as input:
            #        self.__structuresinst = pickle.load(input)
                self.__structuresinst = Object.load(self.__workdir+'/structures.pkl')
                atoms = self.__structuresinst.get_structures()
            else:
                raise Exception("Please do setup first!")
                
        else:
            atoms = self.__structuresinst.get_structures()
        
        for stype, eta in atoms:
            self.__status[(stype,eta)] = {}
            try:
                et.parse(atoms[(stype,eta)].path+'/vasprun.xml')
                self.__status[(stype,eta)]['status'] = 'finished'
                atoms[(stype,eta)].status = 'finished'
            except:
                self.__status[(stype,eta)]['status'] = '--------'
                atoms[(stype,eta)].status = '--------'
            self.__status[(stype,eta)]['path'] = atoms[(stype,eta)].path
        
        return self.__status, self.__structuresinst
    
    def get_filestructure(self):
        return self.__filestructure
    
    