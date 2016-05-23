import os
import time
import threading
import pickle
import json
import Queue as queue
import subprocess

class threads(object):
    def __init__(self, structfile):
        self.__structfile = structfile
        
        
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
        
        self.__vscjobs_dic = dic['cluster']
        
        self.__currpath = None
        
        
        self.__logstr = ""
        
        return
    
    def submit_kstr(self):
        
        self.__logstr += 'copying run_kstr to calculation directory: {0} \n'.format(self.__currpath)
        ## Copy queuing script to calc directory:
        proc = subprocess.Popen(['cp run_kstr {0}'.format(self.__currpath)], shell=True)
        proc.wait()
        workdir = os.getcwd()
        os.chdir(self.__currpath)
        self.__logstr += 'Submitting batch job: {0}/run_kstr \n'.format(self.__currpath)
        proc = subprocess.Popen(['sbatch run_kstr'], shell=True)
        proc.wait()
        os.chdir(workdir)
        return
    
    def submit_shape(self):
        ## Copy queuing script to calc directory:
        self.__logstr += 'copying run_shape to calculation directory: {0} \n'.format(self.__currpath)
        proc = subprocess.Popen(['cp run_shape {0}'.format(self.__currpath)], shell=True)
        proc.wait()
        workdir = os.getcwd()
        os.chdir(self.__currpath)
        self.__logstr += 'Submitting batch job: {0}/run_kstr \n'.format(self.__currpath)
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
        self.__logstr += 'Submitting batch job: {0}/run_kgrn \n'.format(self.__currpath)
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
        self.__logstr += 'Submitting batch job: {0}/run_kfcd \n'.format(self.__currpath)
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
                    
                    self.__logstr += '#####################\n{0} FINISHED \n#####################\n'.format(path)
                    
                    self.__q.task_done()
                    Finished=True
                f.close()
        
        return True
    
    def kstr(self):
        self.submit_kstr()
    
    def start_jobs(self):
        f = open(self.__structfile)
        structobj = pickle.load(f)
        f.close()
        
        
        
        paths=[]
        for struct in structobj.get_structures().values():
            paths.append(struct.path.split('/')[-2]+'/'+struct.path.split('/')[-1]+'/')
        ######################################################
        
        
            
        self.__q=queue.Queue() 
        t=[]   
        
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__kstrpath)
            
            if os.path.exists('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)) and os.stat('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)).st_size != 0: raise SystemExit('Existing calculations in {0}/prn/{1}. Please clean up first!'.format(self.__currpath, self.__kstrname))
            
            self.submit_kstr()
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__kstrname)))
            
        for task in t: 
            task.deamon=True
            task.start()
        self.__q.join()
        print 'KSTR FINISHED --> Starting SHAPE calculations. \n\n'
        
        #### Wait until KSTR has finished.... start SHAPE....
        
        self.__q=queue.Queue() 
        t=[]   
        print paths
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__shapepath)
            
            if os.path.exists('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)) and os.stat('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)).st_size != 0: raise SystemExit('Existing calculations in {0}/prn/{1}. Please clean up first!'.format(self.__currpath, self.__kstrname))
            
            self.submit_shape()
            
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__shapename)))
        for task in t: 
            task.deamon=True
            task.start()
        self.__q.join()
        print 'SHAPE FINISHED --> Starting KGRN calculations. \n\n'
        
        #### Wait until SHAPE has finished.... start KGRN....
        
        self.__q=queue.Queue() 
        t=[]   
        print paths
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__kgrnpath)
            
            if os.path.exists('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)) and os.stat('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)).st_size != 0: raise SystemExit('Existing calculations in {0}/prn/{1}. Please clean up first!'.format(self.__currpath, self.__kstrname))
            
            self.submit_kgrn()
            
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__kgrnname)))
        for task in t: 
            task.deamon=True
            task.start()
        self.__q.join()
        print 'KGRN FINISHED --> Starting KFCD calculations. \n\n'
        
        #### Wait until KGRN has finished.... start KFCD....
        
        self.__q=queue.Queue() 
        t=[]   
        print paths
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__kfcdpath)
            
            if os.path.exists('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)) and os.stat('{0}/prn/{1}'.format(self.__currpath, self.__kstrname)).st_size != 0: raise SystemExit('Existing calculations in {0}/prn/{1}. Please clean up first!'.format(self.__currpath, self.__kstrname))
            
            self.submit_kfcd()
            
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__kfcdname)))
        for task in t: 
            task.deamon=True
            task.start()
        
        
        try:
            self.__q.join()
            print "ALL CALCULATIONS FINISHED."
            
            
        except (KeyboardInterrupt, SystemExit):
            self.__q.put('exit')
            self.__q.join()
            print "Manually terminated!"
        finally:
            with open("log.out",'w') as f:
                f.write(self.__logstr)
        
        
        
#t = threading.Thread(target=subproc, args=(a,))
#t.start()


