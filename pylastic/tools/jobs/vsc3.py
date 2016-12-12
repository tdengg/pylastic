"""
Starting and submitting scripts on vsc3.
"""

import os
import sys
import glob
import time
import threading
import pickle
import json
import Queue as queue
import subprocess



class threads(object):
    def __init__(self, structfile):
        self.__structfile = structfile
        self.__starttime = time.time()
        
        f=open('setup.json')
        dic = json.load(f)
        f.close()
        
        self.__kstrpath = dic['emto']['pnames']['kstr']
        self.__kstrname = '%s.prn'%dic['emto']['jobnames']['structure']
        self.__shapepath = dic['emto']['pnames']['shape']
        self.__shapename = '%s.prn'%dic['emto']['jobnames']['structure']
        self.__kgrnpath = dic['emto']['pnames']['kgrn']
        self.__kgrnname = '%s.prn'%dic['emto']['jobnames']['system']
        self.__kfcdpath = dic['emto']['pnames']['kfcd']
        self.__kfcdname = '%s.prn'%dic['emto']['jobnames']['system']
        
        self.__currpath = None
        return
    
    def submit_kstr(self):
        
        #self.__flog.write('copying run_kstr to calculation directory: {0} \n'.format(self.__currpath))
        ## Copy queuing script to calc directory:
        proc = subprocess.Popen(['cp run_kstr {0}'.format(self.__currpath)], shell=True)
        proc.wait()
        
        workdir = os.getcwd()
        
        os.chdir(self.__currpath)
        self.__flog.write('\t Submitting batch job: {0}/run_kstr \n'.format(self.__currpath))
        proc = subprocess.Popen(['sbatch run_kstr'], shell=True)
        proc.wait()
        
        os.chdir(workdir)
        return
    
    def submit_shape(self):
        ## Copy queuing script to calc directory:
        #self.__flog.write('copying run_shape to calculation directory: {0} \n'.format(self.__currpath))
        proc = subprocess.Popen(['cp run_shape {0}'.format(self.__currpath)], shell=True)
        proc.wait()
        
        workdir = os.getcwd()
        
        os.chdir(self.__currpath)
        self.__flog.write('\t Submitting batch job: {0}/run_kstr \n'.format(self.__currpath))
        proc = subprocess.Popen(['sbatch run_shape'], shell=True)
        proc.wait()
        
        os.chdir(workdir)
        return
    
    def submit_kgrn(self):
        ## Copy queuing script to calc directory:
        
        proc = subprocess.Popen(['cp run_kgrn {0}'.format(self.__currpath)], shell=True)
        proc.wait()
        
        #time.sleep(0.5)
        
        workdir = os.getcwd()
        
        os.chdir(self.__currpath)
        self.__flog.write('\t Submitting batch job: {0}/run_kgrn \n'.format(self.__currpath))
        proc = subprocess.Popen(['sbatch run_kgrn'], shell=True)
        proc.wait()
        
        os.chdir(workdir)
        return
    
    def submit_kfcd(self):
        ## Copy queuing script to calc directory:
        proc = subprocess.Popen(['cp run_kfcd {0}'.format(self.__currpath)], shell=True)
        proc.wait()
        
        workdir = os.getcwd()
        
        os.chdir(self.__currpath)
        self.__flog.write('\t Submitting batch job: {0}/run_kfcd \n'.format(self.__currpath))
        proc = subprocess.Popen(['sbatch run_kfcd'], shell=True)
        proc.wait()
        
        os.chdir(workdir)
        return
    
    def checkstatus(self,path,fname):

        Finished = False
        
        self.__q.put(path)
        
        while not Finished:
            time.sleep(5)
            if os.path.exists('{0}/prn/{1}'.format(path,fname)) and os.stat('{0}/prn/{1}'.format(path,fname)).st_size != 0:
                f=open('{0}/prn/{1}'.format(path,fname))
                lastline = f.readlines()[-1].split()
                if 'Finished' in lastline or 'Volumes:' in lastline: # or 'Volumes:'
                    stime=(time.time()-self.__starttime)
                    M,S=divmod(stime,60)
                    H,M=divmod(M,60)
                    self.__flog.write('#####################\n{0} FINISHED \n Time: {1:02d}:{2:02d}:{3:02d} \n--------------------\n'.format(path, int(H), int(M), int(S)))
                    self.__flog.flush()
                    self.__q.task_done()
                    Finished=True
                for n in range(len(f.readlines())): 
                    if 'Not converged' in f.readlines()[n].split():
                        self.__flog.write('KGRN calculation for {0} NOT CONVERGED..... STOPPING!'.format(path))
                        self.__flog.close()
                        sys.exit('KGRN calculation for {0} NOT CONVERGED..... STOPPING!!'.format(path))
                
                f.close()
            slurmout = glob.glob('{0}/slurm*'.format(path))
            if len(slurmout)!=0:
                with open(slurmout[-1]) as f1:
                    slurmlines = f1.readlines()
                for line in slurmlines:
                    if 'error' in line.split(): 
                        self.__flog.write('ERROR in {0}/prn/{1}!!!! Check slurm output!'.format(path,fname))
                        self.__flog.close()
                        raise SystemExit('ERROR in {0}/prn/{1}!!!! Check slurm output!'.format(path,fname))
        return True
    
    def kstr(self):
        self.submit_kstr()
    
    def start_jobs(self):
        f = open(self.__structfile)
        structobj = pickle.load(f)
        f.close()
        
        self.__flog = open('log.pyl','w')
        
        
        
        paths=[]
        for struct in structobj.get_structures().values():
            paths.append(struct.path.split('/')[-2]+'/'+struct.path.split('/')[-1]+'/')
        ######################################################
        
        
            
        self.__q=queue.Queue() 
        t=[]   
        
        self.__flog.write('Starting KSTR calculations..... \n')
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__kstrpath)
            
            if os.path.exists('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)) and os.stat('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)).st_size != 0: 
                self.__flog.write('Existing calculations in {0}/prn/{1}. Please clean up first!'.format(self.__currpath, self.__kstrname))
                self.__flog.close()
                raise SystemExit('Existing calculations in {0}/prn/{1}. Please clean up first!'.format(self.__currpath, self.__kstrname))
            
            self.submit_kstr()
            self.__flog.flush()
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__kstrname)))
            
        for task in t: 
            task.deamon=True
            task.start()
        self.__q.join()
        self.__flog.write('KSTR FINISHED --> Starting SHAPE calculations. \n\n')
        
        #### Wait until KSTR has finished.... start SHAPE....
        
        self.__q=queue.Queue() 
        t=[]   
        
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__shapepath)
            
            if os.path.exists('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)) and os.stat('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)).st_size != 0: 
                self.__flog.write('Existing calculations in {0}/prn/{1}. Please clean up first!'.format(self.__currpath, self.__kstrname))
                self.__flog.close()
                raise SystemExit('Existing calculations in {0}/prn/{1}. Please clean up first!'.format(self.__currpath, self.__kstrname))
            
            self.submit_shape()
            self.__flog.flush()
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__shapename)))
            
        for task in t: 
            task.deamon=True
            task.start()
        self.__q.join()
        self.__flog.write('SHAPE FINISHED --> Starting KGRN calculations. \n\n')
        
        #### Wait until SHAPE has finished.... start KGRN....
        
        self.__q=queue.Queue() 
        t=[]   
        
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__kgrnpath)
            
            if os.path.exists('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)) and os.stat('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)).st_size != 0: 
                self.__flog.write('Existing calculations in {0}/prn/{1}. Please clean up first!'.format(self.__currpath, self.__kstrname))
                self.__flog.close()
                raise SystemExit('Existing calculations in {0}/prn/{1}. Please clean up first!'.format(self.__currpath, self.__kstrname))
            
            self.submit_kgrn()
            self.__flog.flush()
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__kgrnname)))
        for task in t: 
            task.deamon=True
            task.start()
        self.__q.join()
        self.__flog.write('KGRN FINISHED --> Starting KFCD calculations. \n\n')
        
        #### Wait until KGRN has finished.... start KFCD....
        
        self.__q=queue.Queue() 
        t=[]   
        
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__kfcdpath)
            
            if os.path.exists('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)) and os.stat('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)).st_size != 0: 
                self.__flog.write('Existing calculations in {0}/prn/{1}. Please clean up first!'.format(self.__currpath, self.__kstrname))
                self.__flog.close()
                raise SystemExit('Existing calculations in {0}/prn/{1}. Please clean up first!'.format(self.__currpath, self.__kstrname))
            
            self.submit_kfcd()
            self.__flog.flush()
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__kfcdname)))
            
        for task in t: 
            task.deamon=True
            task.start()
        
        
        try:
            self.__q.join()
            print "ALL CALCULATIONS FINISHED."
            stime=(time.time()-self.__starttime)
            M,S=divmod(stime,60)
            H,M=divmod(M,60)
            self.__flog.write('\n ALL CALCULATIONS FINISHED! \n Total time: {0:02d}:{1:02d}:{2:02d}'.format(int(H), int(M), int(S)))
            
            
        except (KeyboardInterrupt, SystemExit):
            self.__q.put('exit')
            self.__q.join()
            print "Manually terminated!"
        finally:
            self.__flog.close()
            
class vasp(object):
    def __init__(self, structfile):
        self.__structfile = structfile
        self.__starttime = time.time()
        
        f=open('setup.json')
        dic = json.load(f)
        f.close()
        
        
        self.__path = dic['vasp']['pnames']
        self.__name = 'OUTCAR'
        
        self.__currpath = None
        return
    
    
    def submit_vasp(self):
        ## Copy queuing script to calc directory:
        proc = subprocess.Popen(['cp run_vasp {0}'.format(self.__currpath)], shell=True)
        proc.wait()
        
        workdir = os.getcwd()
        
        os.chdir(self.__currpath)
        self.__flog.write('\t Submitting batch job: {0} \n'.format(self.__currpath))
        proc = subprocess.Popen(['sbatch run_vasp'], shell=True)
        proc.wait()
        
        os.chdir(workdir)
        return
    
    def checkstatus(self,path,fname):

        Finished = False
        
        self.__q.put(path)
        
        while not Finished:
            time.sleep(5)
            if os.path.exists('{0}/{1}'.format(path,fname)) and os.stat('{0}/{1}'.format(path,fname)).st_size != 0:
                f=open('{0}/{1}'.format(path,fname))
                lines = f.readlines()
                for line in lines:
                    if 'General timing and accounting informations for this job:' in line: # or 'Volumes:'
                        stime=(time.time()-self.__starttime)
                        M,S=divmod(stime,60)
                        H,M=divmod(M,60)
                        self.__flog.write('#####################\n{0} FINISHED \n Time: {1:02d}:{2:02d}:{3:02d} \n--------------------\n'.format(path, int(H), int(M), int(S)))
                        self.__flog.flush()
                        self.__q.task_done()
                        Finished=True
                
                f.close()
            slurmout = glob.glob('{0}/slurm*'.format(path))
            if len(slurmout)!=0:
                with open(slurmout[-1]) as f1:
                    slurmlines = f1.readlines()
                for line in slurmlines:
                    if 'error' in line.split(): 
                        self.__flog.write('ERROR in {0}. Check slurm output!'.format(path,fname))
                        #self.__flog.close()
                        #raise SystemExit('ERROR in {0}!!!! Check slurm output!'.format(path,fname))
        return True
    

    def start_jobs(self):
        f = open(self.__structfile)
        structobj = pickle.load(f)
        f.close()
        
        self.__flog = open('log.pyl','w')
        
        
        
        paths=[]
        for struct in structobj.get_structures().values():
            paths.append(struct.path.split('/')[-2]+'/'+struct.path.split('/')[-1]+'/')
        ######################################################
        
        
            
        self.__q=queue.Queue() 
        t=[]   
        
        self.__flog.write('Starting calculations..... \n')
        for path in paths:
            self.__currpath = '{0}'.format(path)
            
            if os.path.exists('{0}/{1}'.format(self.__currpath,'vasprun.xml')) and os.stat('{0}/{1}'.format(self.__currpath,'vasprun.xml')).st_size != 0: 
                self.__flog.write('Existing calculations in {0}. Please clean up first!'.format(self.__currpath))
                #self.__flog.close()
                raise SystemExit('Existing calculations in {0}. Please clean up first!'.format(self.__currpath))
            
            self.submit_vasp()
            self.__flog.flush()
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__name)))
            
        for task in t: 
            task.deamon=True
            task.start()
        self.__q.join()
        self.__flog.write('FINISHED. \n\n')
        
        
        
        try:
            self.__q.join()
            print "ALL CALCULATIONS FINISHED."
            stime=(time.time()-self.__starttime)
            M,S=divmod(stime,60)
            H,M=divmod(M,60)
            self.__flog.write('\n ALL CALCULATIONS FINISHED! \n Total time: {0:02d}:{1:02d}:{2:02d}'.format(int(H), int(M), int(S)))
            
            
        except (KeyboardInterrupt, SystemExit):
            self.__q.put('exit')
            self.__q.join()
            print "Manually terminated!"
        finally:
            self.__flog.close()   
        
        
#t = threading.Thread(target=subproc, args=(a,))
#t.start()


