import os
import datetime as dt
import shutil
import errno

months = ['', 'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def createdirs (first_date, tlength, runlength, basedir = "/ocean/sstevens/ariane/"):
    """
    Creates reposities where results and basic infomartion will be stored
    """
    
    runlength = runlength - 1
    
    last_date = first_date + dt.timedelta(hours = 24*runlength)
  

    fdate = '{:%Y%m%d}'.format(first_date) + "_" + '{:%Y%m%d}'.format(last_date) + "_" + "{}d".format(tlength)


    arianedir = basedir + "arianefiles/" + fdate
    resultsdir = basedir + "results/" + fdate
    
    for directory in (arianedir, resultsdir):
    
        if os.path.exists(directory) and directory != '/ocean/sstevens/ariane':
            print ("will delete old version at {}".format(directory))
            shutil.rmtree(directory) #CAUTION HERE!

        os.makedirs(directory)

    
    return {'arianedir': arianedir, 'resultsdir': resultsdir}
    
