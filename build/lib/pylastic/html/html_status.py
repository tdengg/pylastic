import pickle
import json

from lxml import html

class status(object):
    
    def __init__(self, statusdict):
        
        self.__statusdict = statusdict
        
        