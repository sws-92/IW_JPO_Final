import os
import datetime as dt


months = ['', 'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']



def link(n, date, directory, ds):

    for i in ("U", "V", "W", "T"):
        linkfiles(n = n, date = date, directory = directory, vel = i, ds = ds)
        
    print ("velocities linked in", directory)
    
    


def smon (mon):
    months = ['', 'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    return months[mon]


# In[35]:

def sadd0(j):
    if j < 10:
        return "0"+str(j)
 
    else:
        return str(j)
    

def fd(date):
    dd = {
    "day": sadd0(date.day),
    "month": sadd0(date.month),
    "year": sadd0(date.year),
    "smonth": smon(date.month),
    "year2": sadd0(date.year)[-2:]
     }
    return dd


def linkfiles (vel, n, date, directory, ds):
    
    for i in range(1, n+1):

        f = fd(date)

        #source
        source = (f["day"], f["smonth"], f["month"], f["day"], f["month"], f["day"], vel, f["year2"], f["year"], ds )

        src = "{9}/{0}{1}{7}/SalishSea_1h_{8}{2}{3}_{8}{4}{5}_grid_{6}.nc".format(*source)
        
        
        #destination
        dst = directory + "/SalishSea_{:03d}_grid_{}.nc".format(i, vel)

        if not isinstance(i,int):
            raise ValueError('FAIL!!!!!!!!!!!!!')

        os.symlink(src, dst)

        date = date + dt.timedelta(hours = 24)


