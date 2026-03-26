# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 12:00:14 2026
Last Update : Fri Mar 6 2026
@author: goirand-lopez
"""
import os
from os.path import join


# !!! Before using the code you have to select the right folders if you don't want in the same folder !
cellpath = os.getcwd() #  Path of the folder where your neuron models are
scriptpath = os.getcwd() #  Path of the folder where your script is
resultsfolderpathname = cellpath # Path of the folder where you want to store the results

PPnspk = 10 # Number of GC activated by the initial stimulation, 5 or 10 in our model
PPWindowsLength = 0 # Duration of the regular stimulation windows [0, 15, 30, 45, 60, 75, 90, 105] ms
PPint = PPWindowsLength/PPnspk # Calculated regular intervals between initial stimulation
nspr = 10 # nspr among [10, 13, 16, 19, 22, 25, 28, 31]
IP = 0 # IP among [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    
ParametersSimulation = dict(
nGClayer = 25, # Number of layer in GCL 5 or 25 in our models
PPnspk = PPnspk , # Number of GC activated by the initial stimulation, 5 or 10 in our model
PPintervals = PPint, # PPint is calculated to have a PPnspk(5 or 10) regular intervals between sitmulation in windows of [0, 15, 30, 45, 60, 75, 90, 105] ms
simvar = True, # Do we simulate the activity of not
nspr = nspr, # nspr among [10, 13, 16, 19, 22, 25, 28, 31]
lightstim = False, # Addition or not of a light stimulation during the simulation after initial stimulation
GCstim = 1, # Number of GCs stimulated by each PP,
IP = 0,# IP among [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
randseed = False, # Activation of randomized the seed used to generate random events
saveopt = True # Option choice to save or not results
)
    
with open(os.path.join(cellpath, 'ParametersSimulationSingle'+'.txt'), "w") as f:
    for key,val in ParametersSimulation.items() :
        f.write('{}\t{}\t{}\t'.format(key,val,type(val)))
        f.write('\n')

exec(open(join(scriptpath,'Network500PPs_1IPfor1window_1PP1GC_singleprocess.py')).read())