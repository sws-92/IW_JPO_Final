import datetime as dt
from shutil import copyfile
import inspect
from create_dirs import createdirs
from create_namelist import namelist
from create_inpos import longinitp
from link_vels import link as link
from run_ariane import run_ariane

def pre_run (first_day, last_day, trajectory_length, src, ds):

    """
    Prepares the run and starts the run (calls start_run)

    :arg first_day: initial day of the run
    :type to: datetime object

    :arg trajectory_length: trajectory length of each particle in days
    :type tf: positive integer
    
    :arg ds: velocties data source (eg, "nowcast", "hindcast")
    :type tf: string

    :returns: this function does not have a return value
    """

    run_length = (last_day - first_day).days + 1

    #starts the runs
    start_run(first_day = first_day, 
              last_day = last_day, 
              run_length = run_length, 
              trajectory_length = trajectory_length,
              src_code = src,
              ds = ds)
    
def start_run(
        first_day, 
        last_day, 
        trajectory_length, 
        run_length,
        src_code,
        ds):
    '''
    starts the run and collects the results
    '''
    
    #creates directories
    dirs = createdirs(first_day, trajectory_length, run_length)
    arianedir = dirs["arianedir"]
    resultsdir = dirs["resultsdir"]
    

    #creates namelist
    namelist(arianedir, trajectory_length + run_length, trajectory_length)
    
    
    #creates initial_positions.txt
    longinitp(arianedir, run_length )
    
    
    #links velocties (.nc files, U, V, W as default)
    link(trajectory_length + run_length , first_day, arianedir, ds)


    #copies basic files
    src = arianedir + "/namelist"
    dst = resultsdir + "/namelist"
    copyfile(src, dst)

    src = arianedir + "/initial_positions.txt"
    dst = resultsdir + "/initial_positions.txt"
    copyfile(src, dst)
    
    #copies code that originated the run to the results directory
    src = src_code
    dst = resultsdir + "/source.py"
    copyfile(src, dst)
    
    #runs ariane
    time_before = dt.datetime.now()
    error, log, final_message = run_ariane(arianedir, resultsdir)
    time_after = dt.datetime.now()

    #copies results to the results directory
    src = arianedir + "/traj.txt"
    dst = resultsdir + "/traj.txt"
    copyfile(src, dst)

    #writes the information about the run and adds it to the respective result directory
    run_info = {
                "first_day" : first_day, 
                "last_day" : last_day, 
                "time_before": time_before.strftime("%Y-%m-%d %H:%M:%S"), 
                "time_after" : time_after.strftime("%Y-%m-%d %H:%M:%S"), 
                "time_delta" : time_after - time_before, 
                "source" : src_code,
                "ds": ds,
                "traj_len": trajectory_length
                }
    
    run_message = (
    "The simulated run starts on {first_day} and finishes on {last_day}.\n"
    "Each trajectory is {traj_len}-day long.\n"
    "This sumulation was performed from {time_before} to {time_after}, taking {time_delta} to complete.\n"
    "The velocities are linked from {ds} data.\n"
    "Please refer to {source} to see the code that originated these results (also available in this directory).\n"
                  ).format (**run_info)
        
    with open(resultsdir + "/run_info", 'w') as file:
        file.write(run_message.format(**run_info)) 

    with open(resultsdir + "/run_log.txt", 'w') as file:
        file.write(log)
