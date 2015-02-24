import pickle

class Object(object):
    def __init__(self):
        self.__object = None
    
    def load(self, f):
        with open(f, 'rb') as input:
            self.__object = pickle.load(input)
        return self.__object
    
    def save(self, f, instance):
        with open(f, 'wb') as output:
            pickle.dump(instance, output, -1)
        return
    
    