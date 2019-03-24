import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
from datetime import datetime

def mode_choose2(file, mode, n, N, rho, T):
    
    current_mode = (file + ' ' +
                    mode + ' ' +
                    str(n) + ' ' +
                    str(N) + ' ' +
                    str(rho) + ' ' +
                    str(T)
                   )

    return(current_mode)

def mode_choose3(file, mode, n, N, rho, T, rCutOff, rMax, updateInterval):
    
    current_mode = (file + ' ' +
                    mode + ' ' +
                    str(n) + ' ' +
                    str(N) + ' ' +
                    str(rho) + ' ' +
                    str(T) + ' ' + 
                    str(rCutOff) + ' ' +
                    str(rMax) + ' ' +
                    str(updateInterval)
                   )

    return(current_mode)

# MAIN

# Running through various rCutOffs, rMaxs, updateIntervals
i_min = 2.5
i_max = 3.5

j_min = 0.5
j_max = 2.5

k_min = 5
k_max = 15

i_s = (i_max - i_min) * 10
j_s = (j_max - j_min) * 10
k_s = k_max - k_min
number_of_plots = k_s * (j_s + 1) * (i_s + 1)

# Simulated steps
n = 400

# Width of progressbar on screen
progressbar_width = 20

# Count all steps
counting = 0

# Measur runtime of full dump
starttime = datetime.now()
with open('runtimes.dat', 'w+') as dataFile:

    current_mode = mode_choose2(file='..\Release\md2.exe', mode='bounded', n=n, N=32, rho=0.95, T=1.0)
    
    # Run init. velocity-controled simulation with current conditions and measure its runtime
    current_step_starttime = datetime.now()
    os.system(current_mode)
    dataFile.write('{0}\t{1}\t{2}\t{3}\n'.format(0, 0, 0, (datetime.now() - current_step_starttime).total_seconds()))

    for i in np.linspace(i_min, i_max, i_s+1):
        for j in np.linspace(i + j_min, i + j_max, j_s+1):
            for k in range(k_min, k_max):
                
                current_mode = mode_choose3(file='..\Release\md3.exe', mode='bounded', n=n, N=32, rho=0.95, T=1.0, rCutOff=i, rMax=j, updateInterval=k)
                
                # Run FULL simulation with current conditions and measure its runtime
                current_step_starttime = datetime.now()
                os.system(current_mode)
                dataFile.write('{0}\t{1}\t{2}\t{3}\n'.format(i, j, k, (datetime.now() - current_step_starttime).total_seconds())) 

                # Count runs
                counting += 1

                sys.stdout.write("\rProgress: [{0}{1}] {2:.3f}%\tElapsed: {3}\tRemaining: {4}\t{5}".format(u'\u2588' * int(counting/number_of_plots * progressbar_width),
                                                                                                           u'\u2591' * (progressbar_width - int(counting/number_of_plots * progressbar_width)),
                                                                                                           counting/number_of_plots * 100,
                                                                                                           datetime.now()-starttime,
                                                                                                           (datetime.now()-starttime)/counting * (number_of_plots - counting),
                                                                                                           counting))
                sys.stdout.flush()