import numpy as np
import __main__

class CQP(object):
    def __init__(self, ibrav='bcc'):
        self.ibrav = ibrav
        
        self.struc = {'bcc':{
            'b_1':2.*np.pi*np.array([1.,1.,0.]), 'b_2':2.*np.pi*np.array([0.,1.,1.]),'b_3':2.*np.pi*np.array([1.,0.,1.]),
            'a_1':1./2.*np.array([1.,1.,-1.]), 'a_2':1./2.*np.array([-1.,1.,1.]),'a_3':1./2.*np.array([1.,-1.,1.]),
            'ac_1':np.array([1.,0.,0.]), 'ac_2':np.array([0.,1.,0.]),'ac_3':np.array([0.,0.,1.])
                 }
        }
    def calc_cqp(self):
        k=np.dot(self.struc['bcc']['b_1'],self.struc['bcc']['ac_1'])/np.pi
        
        k_uvw = lambda x: x[0]*self.struc['bcc']['b_1'] + x[1]*self.struc['bcc']['b_2'] + x[2]*self.struc['bcc']['b_3']
        
        #points = np.linspace(0,1,1001)
        points=[1]
        
        for N in [2,3,4,5,6,7,8]:
            for p in points:
                
                #uvw = [-0.5*p, 0.5*p, 0.5*p]
                #uvw = [-0.5*(1-p*0.75), 0.5*(1-p*0.25), 0.5*(1-p*0.25)]
                uvw = [0.,0.5,0.5]
                n_xyz = [N,N,N]
                
                for i in range(3):
                    a_n = n_xyz[i]*self.struc['bcc']['a_%s'%(i+1)]
        
                    adotk = (np.dot(a_n,k_uvw(uvw)))

                    if adotk/(2.*np.pi)%1==0 and adotk/(2.*np.pi)>0.: 
                        print adotk/(2.*np.pi), uvw, N
                
        
        
        return k
        
if __name__ == "__main__":
    print 'HALLO'
    c = CQP()
    c.calc_cqp()