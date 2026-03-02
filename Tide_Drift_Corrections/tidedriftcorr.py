#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: bronsoncvijanovich
"""

import numpy as np
import tamura
import pandas as pd
from datetime import datetime 

# %% Data import and extraction

fullsheet = pd.read_csv('TIMP_2024_allgravitydata_01282025.csv') # import the CG6 csv data

stations = fullsheet['Station']
dates = fullsheet['Date'] # define data columns from csv
times = fullsheet['Time']
lat = fullsheet['Lat_rtk']
lon = fullsheet['Lon_rtk']
elev = fullsheet['Elev_rtk']
g = fullsheet['CorrGrav'] # this is called CorrGrav because it is corrected for the CG6's sensor temp and tilt - we kept the drift and tide corrections off because we can calculate them moe accurately here
g_abs_fasb = 979770.20147 # +/- 0.00184 this is the absolute gravity measurement taken by the FG5 at FASB at UofU in SLC in mGal


dateandtime = np.empty([len(dates),6]) # arrray with columns: year, month, day, hour, minute, second for all grav measurements
dt = np.empty(len(dates)).astype(datetime)
for i in range(0,len(dates)):        # this for loop puts dates and times from the CG6 into an array
    date_sep = dates[i].split('-')
    dateandtime[i,0] = date_sep[0] # year column
    dateandtime[i,1] = date_sep[1] # month column
    dateandtime[i,2] = date_sep[2] # day column
    time_sep = times[i].split(':')
    dateandtime[i,3] = time_sep[0] # hour column
    dateandtime[i,4] = time_sep[1] # minute column
    dateandtime[i,5] = time_sep[2] # second column
    dt[i] = datetime(int(date_sep[0]),int(date_sep[1]),int(date_sep[2]),int(time_sep[0]),int(time_sep[1]),int(time_sep[2])) # column vector of the datetimes of all meausrements to allow for elapsed time calculations in drift correction section
    
dateandtime = dateandtime.astype(int) # integerize date and time array to allow datetime function to work below

# %% Tide Correction  
tl = 0
ot = 0    # conversion from local time
tide_corr = np.zeros(len(dateandtime))

for i in range (0,len(dateandtime)): # calculate the tide correction for each data point
    tide_corr[i] = (tamura.tide(dateandtime[i,0], dateandtime[i,1], dateandtime[i,2], dateandtime[i,3], dateandtime[i,4], dateandtime[i,5], lon[i], lat[i], elev[i], tl, ot))/1000


g = g + tide_corr # apply tide correction to gravity data

# %% At-Timp Drift Correction

drift_rows = [[3,14],[22,36],[37,50],[51,61],[67,78],[80,92],[93,104],[105,116],[117,128],[129,141],[142,153],[163,174],[174,187],[188,201],[202,213],[213,220],[226,232],[232,243],[243,256],[256,260],[266,277],[277,285],[286,296],[296,307]] # each pair of numbers is the row number of two drift measurements that bracket a group of gravity measurements on the glacier
tgbs_benchmark = g[drift_rows[0][0]]

drift_rates = np.zeros(len(drift_rows)) # initialize array to store drift rates
drift_corrections = np.zeros(309) # initialize array to store drift corrections

for i in range(0,len(drift_rows)):
    start = drift_rows[i][0] # row of starting drift measurement
    end = drift_rows[i][1] # row of ending drift measurement
    data_rows = np.arange(start+1,end) # the row numbers of the gravity measurements for this drift block
    g_start = g[start] # starting drift measurement (taken at base station before going onto the rock glacier)
    g_end = g[end] # ending drift measurement (taken at base station upon return from rock glacier)
    start_time = dt[start] # time of starting drift measurement
    end_time = dt[end] # time of ending drift measurement
    time_diff = (end_time-start_time).total_seconds() # length of time between the two measurements in seconds
    drift_rate = (g_end-g_start)/time_diff # drift rate in mGal/s
    drift_rates[i] = drift_rate # store drift rates in array
    g[end] = g[end] - drift_rate*time_diff # correct ending drift measurement
    timp_shift = g[drift_rows[i][0]] - tgbs_benchmark # difference in g between drift blocks
    g[start] = g[start] - timp_shift # applying dc shift to starting drift measurement
    g[end] = g[end] - timp_shift # applying dc shift to ending drift measurement
    drift_corrections[drift_rows[i][1]] = drift_rate*time_diff*-1 # store drift cirrections
    for j in data_rows: # iterate through this block of measurements and correct them for drift
        elapsed = (dt[j]-start_time).total_seconds() # elapsed time between gravity measurement and starting drift measurement
        drift_corr = drift_rate*elapsed # calculate drift correction for the gravity measurement
        drift_corrections[j] = drift_corr*-1 # store drift corrections
        g[j] = g[j] - drift_corr - timp_shift # apply drift correction and dc shift to gravity data

# %% FCalculate offsets between absolute gravity base station at FASB and Timp base station (TGBS)
fasb_rows = [[0,16],[17,63],[65,156],[159,222],[223,262],[263,309]] # row indices of FASB measurements
tgbs_end_rows = [[3,14],[22,61],[67,153],[163,220],[226,260],[266,307]] # row indices of survey start/end measurements at the rock glacier base station
survey_num = [1,2,3,4,5,6] # survey number (6 surveys in Fall 2024)
f_t_updiff = np.zeros(6) # initialize arrays for difference in gravity between FASB and TGBS
f_t_downdiff = np.zeros(6)


for i in range(0,len(fasb_rows)):
    f_t_updiff[i] = g[fasb_rows[i][0]]-g[tgbs_end_rows[i][0]] # calculate diff in grav from FASB to TGBS
    f_t_downdiff[i] = g[fasb_rows[i][1]]-g[tgbs_end_rows[i][1]] # calculate diff in grav from TGBS to FASB

ftuavg = np.average(f_t_updiff)*np.ones(6) #find averages
ftdavg = np.average(f_t_downdiff)*np.ones(6)
ft_avg = np.average([ftuavg,ftdavg])


ftustd = np.std(f_t_updiff)*np.ones(6) # calculate standard deviations to check consistency
ftdstd = np.std(f_t_downdiff)*np.ones(6)

# %% Absolute Gravity Tie-in

for i in range(0,len(drift_rows)):
    g[drift_rows[i][0]] = g_abs_fasb - ft_avg
    g[drift_rows[i][1]] = g_abs_fasb - ft_avg
    
    start = drift_rows[i][0] # row of starting drift measurement
    end = drift_rows[i][1] # row of ending drift measurement
    data_rows = np.arange(start+1,end) # the row numbers of the gravity measurements for this drift block
    for j in data_rows:
        g_diff = g[j] - tgbs_benchmark
        g[j] = g[start] + g_diff

g.to_csv('g.csv')
