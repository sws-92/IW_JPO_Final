#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 12:01:56 2017

@author: gsgarbi
"""

import datetime as dt
from dateutil.relativedelta import relativedelta


from start_run import pre_run

'''
Runs ariane up to the most recent avaliable data.
This codes assumes that data up to [(present day) - 2 days] 
is available on nowcast, hindcast or other.

'''

if __name__ == "__main__":
    
    for m in range (8,9):
    
        DATA_LIMIT = dt.date.today() - relativedelta(days = 3)
    
    
        first_day = dt.date(
                            year = 2016,
                            month = m,
                            day = 1
                            )
                 
        
        last_day = first_day + relativedelta(months = 1) - relativedelta(days = 1) #up to the last day of that month
        
    
        
        trajectory_length = min((DATA_LIMIT - last_day).days, 1)
        
        


    
        #print (first_day, last_day, DATA_LIMIT, trajectory_length)
    
    
        
        pre_run (
                     first_day = first_day, 
                     last_day = last_day, 
                     trajectory_length = trajectory_length, 
                     src = '/ocean/sstevens/ariane/1create_repos_run_ariane/customized_run.py',
                     ds = 'hindcast'
                )        
