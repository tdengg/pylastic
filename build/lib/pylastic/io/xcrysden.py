import os
import re
import sys
import numpy as np

class XCD(object):
    """ Module for generating xcrysden input files from atoms object.
    """
    def __init__(self, verbouse=False):
        self.__verbouse=verbouse
        self.__basevect=[]
        self.__lattvect_pos=None
        self.__natoms = 0
        self.__pos={}
        
    def trans_csystem(self, in_mat, out_mat, in_vec):
        """Transformation of coordinate system.
        
        :param 3x3-matrix in_mat: Old coordinate frame.
        
        :param 3x3-matrix out_mat: New coordinate frame.
        
        :param 3-vector in_vec: Vector to transform.
        
        .. math::
            v_{out}= M_{in}^T M_{out}^T v_{in}
        """
        c_vec = np.dot(np.dot(np.linalg.inv(np.transpose(in_mat)),np.transpose(out_mat)), in_vec)
        #c_vec = np.dot(np.dot(np.linalg.inv(in_mat),out_mat), in_vec)
        return c_vec
        
    def write_xcrysden(self, pos, fname):
        path = fname.rstrip(os.path.basename(fname))
        fname = os.path.basename(fname)
        
        os.system('mkdir %s/xcrysden'%(path))
        path=path+'/xcrysden/'
        
        xsf=open((path+fname),'w')
        
        self.__xsf = []
        self.__xsf.append('CRYSTAL \n')
        self.__xsf.append('PRIMVEC \n')
        
        self.__pos=pos
        BS1 = self.__pos['vlatt_1']
        BS2 = self.__pos['vlatt_2']
        BS3 = self.__pos['vlatt_3']
        
        a=np.sqrt(BS1[0]**2.+BS1[1]**2.+BS1[2]**2.)
        b=np.sqrt(BS2[0]**2.+BS2[1]**2.+BS2[2]**2.)
        c=np.sqrt(BS3[0]**2.+BS3[1]**2.+BS3[2]**2.)
        
        scale = self.__pos['scale']
        
        self.__xsf.append('%s %s %s \n'%(BS1[0]*scale, BS1[1]*scale, BS1[2]*scale))
        self.__xsf.append('%s %s %s \n'%(BS2[0]*scale, BS2[1]*scale, BS2[2]*scale))
        self.__xsf.append('%s %s %s \n'%(BS3[0]*scale, BS3[1]*scale, BS3[2]*scale))
        
        B_matrix = np.array([BS1,BS2,BS3])
        C_matrix = np.identity(3)
        
        
        
        QX=[]
        QY=[]
        QZ=[]
        print self.__pos['vbasis']
        C_vec = []
        for i in range(len(self.__pos['natoms'])):
            B_vec = np.array(self.__pos['vbasis']['b_%s'%(i+1)])
            C_vec.append(self.trans_csystem(C_matrix, B_matrix, B_vec)) #Convert direct to cartesian coordinates
        
        for vec in C_vec:
            QX.append(vec[0])
            QY.append(vec[1])
            QZ.append(vec[2])
        
        self.__xsf.append('PRIMCOORD \n')
        self.__xsf.append('%s 1 \n'%len(self.__pos['natoms']))
        for i in range(len(QX)):
        
            self.__xsf.append('%s %s %s %s \n'%(i, QX[i]*scale, QY[i]*scale, QZ[i]*scale))


        
        
        
        
        alpha = np.arccos(np.dot(BS2,BS3)/(b*c))*180./np.pi
        beta =  np.arccos(np.dot(BS3,BS1)/(a*c))*180./np.pi
        gamma = np.arccos(np.dot(BS1,BS2)/(a*b))*180./np.pi
        
        #self.__kstr = self.replace(self.__kstr,"kstr","BSX", BSX)
        #self.__kstr = self.replace(self.__kstr,"kstr","BSY", BSY)
        #self.__kstr = self.replace(self.__kstr,"kstr","BSZ", BSZ)
        #print self.__kstr
        xsf.writelines(self.__xsf)
        xsf.close()
        
        return a
    