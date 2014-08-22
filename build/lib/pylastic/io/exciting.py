import lxml.etree as et
import os

class INPUT(object):
    
    def __init__(self, fname='input.xml'):
        self.__fname = fname
        f = open(self.__fname, 'r')
        self.__doc  = et.parse(f)
        self.__root = self.__doc.getroot()
        
    
    def read_in(self):
        
        #####       Generate spacegroup output       ######
        os.system('xsltproc $ElaSticROOT/exciting2sgroup.xsl '+ self.__fname +' > sgroup.in')
        os.system('sgroup sgroup.in 1>sgroup.out 2>sgroup.err')
        
        
        
        i_dict = {}
        i_dict['name'] = self.__doc.xpath('//title')[0].text
        i_dict['path'] = self.__fname
        i_dict['scale'] = map(float,self.__doc.xpath('/input/structure/crystal/@scale'))[0]
        bv_dump = self.__doc.xpath('//basevect/text()')
        bv = []
        for basevect in bv_dump:
            bv.append(map(float,basevect.split()))
        
        #####  Get spacegroup number from sgroup.out  ######
        SGf   = open('sgroup.out', 'r')
        SGlins= SGf.readlines()
        SGf.close()

        for i in range(len(SGlins)):
            if (SGlins[i].find('Number and name of space group:') >= 0):
                SGN = int(float(SGlins[i].split()[6]))
                SGN_explanation=SGlins[i].strip()
                break
        ####################################################
        i_dict['sgn'] = SGN
        i_dict['vlatt_1'] = bv[0]
        i_dict['vlatt_2'] = bv[1]
        i_dict['vlatt_3'] = bv[2]
        i_dict['natoms'] = len(self.__doc.xpath('//atom'))      #FIX: coherence with other codes
        i_dict["selective"] = None
        i_dict["vbasis"] = {}
        #for i in range(len(p_dict["natoms"])):
        #    p_dict["vbasis"]["species_"+str(i+1)] = []
        #    for j in range(p_dict["natoms"][i]):
        #        p_dict["vbasis"]["species_"+str(i+1)].append(self.lta())
        return i_dict
    
    def write_in(self, i_dict, fileName):
        
        f   = open(fileName, 'w')
        f.write(et.tostring(self.__root, method  ='xml',
                                       pretty_print   =True ,
                                       xml_declaration=True ,
                                       encoding       ='UTF-8'))
        f.close()
        
        return
    
    def read_sgroup(self):
        return
    
    def write_sgroup(self):
        return