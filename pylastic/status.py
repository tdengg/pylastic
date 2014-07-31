import os
import sys
import pickle
import lxml.etree as et

class Check(object):
    def __init__(self):
        self.__starus = None
        self.__structuresinst = None
        self.__workdir = './'
    def check_calc(self):
        self.__status = {}
        self.__paths = {}
        if not self.__structuresinst:
            try:
                with open(self.__workdir+'/structures.pkl', 'rb') as input:
                    atoms = pickle.load(input)
            except:
                raise Exception("Please do setup first!")
                
        else:
            atoms = self.__structuresinst.get_structures()
        
        for stype, eta in atoms:
            self.__status[(stype,eta)] = {}
            try:
                et.parse(atoms[(stype,eta)].path+'/vasprun.xml')
                self.__status[(stype,eta)]['status'] = 'finished'
            except:
                self.__status[(stype,eta)]['status'] = '--------'
            self.__status[(stype,eta)]['path'] = atoms[(stype,eta)].path
        
        return self.__status
    
    def get_filestructure(self):
        return self.__filestructure
    
    