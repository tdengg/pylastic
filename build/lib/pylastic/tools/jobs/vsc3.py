import os
import time
import threading
import pickle
import Queue as queue
import subprocess

class threads(object):
    def __init__(self, structfile):
        self.__structfile = structfile
        self.__kstrpath = 'kstr'
        self.__kstrname = 'NiTi.prn'
        self.__shapepath = 'shape'
        self.__kgrnpath = 'kgrn'
        self.__kfcdpath = 'kfcd'
        self.__currpath = None
        return
    
    def submit_kstr(self):
        print 'sbatch {0}/run_kstr'.format(self.__currpath)
        ## Copy queuing script to calc directory:
        subprocess.Popen(['cp run_kstr {0}'.format(self.__currpath)])
        
        subprocess.Popen(['sbatch {0}/run_kstr'.format(self.__currpath)])
        
        
        return
    
    def submit_shape(self):
        ## Copy queuing script to calc directory:
        subprocess.Popen(['cp run_shape {0}'.format(self.__currpath)])
        
        subprocess.Popen(['sbatch {0}/run_shape'.format(self.__currpath)])
        return
    
    def submit_kgrn(self):
        ## Copy queuing script to calc directory:
        subprocess.Popen(['cp run_kgrn {0}'.format(self.__currpath)])
        
        subprocess.Popen(['sbatch {0}/run_kgrn'.format(self.__currpath)])
        return
    
    def submit_kfcd(self):
        ## Copy queuing script to calc directory:
        subprocess.Popen(['cp run_kfcd {0}'.format(self.__currpath)])
        
        subprocess.Popen(['sbatch {0}/run_kfcd'.format(self.__currpath)])
        return
    
    def checkstatus(self,path,fname):

    
        Finished = False
        
        self.__q.put((path,'done'))
        
        while not Finished:
            time.sleep(5)
            if os.path.exists('{0}/prn/{1}'.format(path,fname)) and os.stat('{0}/prn/{1}'.format(path,fname)).st_size != 0:
                f=open('{0}/prn/{1}'.format(path,fname))
                if 'Finished' or 'Volumes:' in f.readlines()[-1].split(): 
                    Finished=True
                    print '{0} FINISHED'.format(path)
                    
                f.close()
        self.__q.task_done()
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
        
        self.__q=queue.Queue() 
        t=[]   
        print paths
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__kstrpath)
            self.submit_kstr()
            
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__kstrname)))
        for task in t: 
            task.deamon=True
            task.start()
        self.__q.join()
        #### Wait until KSTR has finished.... start SHAPE....
        
        self.__q=queue.Queue() 
        t=[]   
        print paths
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__shapepath)
            self.submit_shape()
            
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__shapename)))
        for task in t: 
            task.deamon=True
            task.start()
        self.__q.join()
        #### Wait until SHAPE has finished.... start KGRN....
        
        self.__q=queue.Queue() 
        t=[]   
        print paths
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__kgrnpath)
            self.submit_kgrn()
            
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__kgrnname)))
        for task in t: 
            task.deamon=True
            task.start()
        self.__q.join()
        print 'FINISHED --> PROCEED'
        
        #### Wait until KGRN has finished.... start KFCD....
        
        self.__q=queue.Queue() 
        t=[]   
        print paths
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__kfcdpath)
            self.submit_kfcd()
            
            t.append(threading.Thread(target=self.checkstatus, args=(self.__currpath, self.__kfcdname)))
        for task in t: 
            task.deamon=True
            task.start()
        self.__q.join()
        
        print "CALCULATIONS FINISHED"
        
#t = threading.Thread(target=subproc, args=(a,))
#t.start()


