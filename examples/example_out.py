import sys
sys.path.append('/home/t.dengg/git/pylastic/pylastic')

from elatoms import Structures, ElAtoms
from postprocess import ECs
from status import Check

import matplotlib.pyplot as plt 


ec = ECs()
ec.set_structures()
ec.set_gsenergy()
ec.set_analytics()

#print ec.get_CVS()
ec.plot_cvs()

ec.set_ec('0.05') 

#print ec.get_ec()

#print(c)
f = ec.plot_2nd()
print ec.get_rms()
plt.show()


