import analyze
import get_DFTdata
import json
import os

class ECs(object):
    def __init__(self):
        self.V0 = None
        
    def set_CVS(self):
        
        getData = get_DFTdata.VASP()
        f=open('info.json')
        dic = json.load(f)
        
        for key in sorted(dic.keys()):
            energy = []
            strain = []
            for key2 in sorted(map(int,dic[key].keys())):
                
                getData.set_outfile('Dst%.2d'%int(key)+'/Dst%.2d'%int(key)+'_%.2d'%int(key2)+'/vasprun.xml')
                getData.set_gsEnergy()
                energy.append(getData.get_gsEnergy())
                strain.append(dic[key][str(key2)]['eta'])
            self.V0 = dic[key][str(key2)]['V0']
            
            ans = analyze.Energy(strain,energy,self.V0)
            ans.set_2nd(6)
            print ans.get_2nd()
            
if __name__ == '__main__':
    ECs().set_CVS()
            
        