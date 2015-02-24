class PrettyMatrix(object):
    def __init__(self):
        
        self.__defMatrix = []
    
    def printListAsMatrix(self, lst, flp):
        prettyMatrix = ''
        i=1
        for element in lst:
            prettyMatrix = prettyMatrix + "{:10d}".format(element)
            if i%3==0: prettyMatrix = prettyMatrix + '\n'
            i+=1
        return prettyMatrix 
    
    def printMatrix(self, matrix, flp):
        
        prettyMatrix = matrix
        return prettyMatrix
    
class FileStructure(object):
    def __init__(self):
        self.filelist = []
        
    def tree(self):
        return
    
    def dicToTree(self, dic):
        temp = []
        string = ''
        for key in sorted(dic.keys()):
            if not key[0] in temp:
                temp.append(key[0])
                string += '%s\n'%str(key[0])
            else:
                string += '\t%s\n'%str(key[1])
                for k in dic[key].keys():
                    string += '\t\t%s\n'%str(dic[key][k])
        
        return string
    
    
    
    