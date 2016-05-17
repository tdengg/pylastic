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
        self.__currpath = None
        
        self.__f_log=open('jobs.log','w')
        
        return
    
    def submit_kstr(self):
        print 'sbatch {0}/run_kstr'.format(self.__currpath)
        print 'cp run_kstr {0}'.format(self.__currpath)
        ## Copy queuing script to calc directory:
        proc = subprocess.Popen(['cp run_kstr {0}'.format(self.__currpath)], shell=True)
        proc.communicate()
        workdir = os.getcwd()
        os.chdir(self.__currpath)
        proc = subprocess.Popen(['sbatch run_kstr'], shell=True)
        proc.communicate()
        os.chdir(workdir)
        return
    
    def submit_shape(self):
        ## Copy queuing script to calc directory:
        proc = subprocess.Popen(['cp run_shape {0}'.format(self.__currpath)], shell=True)
        proc.communicate()
        workdir = os.getcwd()
        os.chdir(self.__currpath)
        proc = subprocess.Popen(['sbatch run_shape'], shell=True)
        proc.communicate()
        os.chdir(workdir)
        return
    
    def submit_kgrn(self):
        ## Copy queuing script to calc directory:
        proc = subprocess.Popen(['cp run_kgrn {0}'.format(self.__currpath)], shell=True)
        proc.communicate()
        workdir = os.getcwd()
        os.chdir(self.__currpath)
        proc = subprocess.Popen(['sbatch run_kgrn'], shell=True)
        proc.communicate()
        os.chdir(workdir)
        return
    
    def submit_kfcd(self):
        ## Copy queuing script to calc directory:
        proc = subprocess.Popen(['cp run_kfcd {0}'.format(self.__currpath)], shell=True)
        proc.communicate()
        workdir = os.getcwd()
        os.chdir(self.__currpath)
        proc = subprocess.Popen(['sbatch run_kfcd'], shell=True)
        proc.communicate()
        os.chdir(workdir)
        return
    
    def checkstatus(self,path,fname):

    
        Finished = False
        
        self.__q.put(path)
        
        while not Finished:
            time.sleep(5)
            if os.path.exists('{0}/prn/{1}'.format(path,fname)) and os.stat('{0}/prn/{1}'.format(path,fname)).st_size != 0:
                f=open('{0}/prn/{1}'.format(path,fname))
                if 'Finished' or 'Volumes:' in f.readlines()[-1].split(): 
                    Finished=True
                    print '{0} FINISHED'.format(path)
                    self.__f_log.write('Finished kstr %s.\n'%path)
                    self.__q.task_done()
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
            
        self.__q=queue.Queue() 
        t=[]   
        
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__kstrpath)
            self.submit_kstr()
            self.__f_log.write('Submitted kstr %s.\n'%path)
            
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
        
        
        try:
            self.__q.join()
            print "CALCULATIONS FINISHED"
            self.__f_log.write("CALCULATIONS FINISHED")
            self.__f_log.close()
        except (KeyboardInterrupt, SystemExit):
            self.__q.put('exit')
            self.__q.join()
            print "Manually terminated!"
            self.__f_log.close()
        
        
#t = threading.Thread(target=subproc, args=(a,))
#t.start()


