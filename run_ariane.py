import subprocess
import os
import signal

def handler(signum, frame):
    raise OSError

def run_ariane(arianedir, resultsdir):
    """This function checks for horizontal and vertical eddies, which are errors that Ariane misses
    It does not catch Ariane internal errors
    It takes as arguments two directory paths and returns a (Boolean, text file, string)"""
    
    os.chdir(arianedir)
    ar = subprocess.Popen('/ocean/sstevens/ariane/arianeCode/ariane-2.2.8-code/bin/ariane', stdout=subprocess.PIPE)
    reader = ar.stdout
    log = '' 

    while ar.poll() is None:

        signal.signal(signal.SIGALRM, handler)
        try:
            signal.alarm(3600)
            # This reader.readline() may hang indefinitely for Horizontal or Vertical eddies  
            w = reader.readline().decode('ascii')
            print (w)
            log = log + w
            signal.alarm(0) #deactivate the alarm
            final_message = w
        except OSError:
            w = reader.readline().decode('ascii')
            log = log + w
            print (w)
            ar.kill()
            return True, log, w
    
    final_message = w
        
    return False, log, final_message
        
