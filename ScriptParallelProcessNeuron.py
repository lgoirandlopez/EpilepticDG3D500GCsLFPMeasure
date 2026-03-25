# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 00:44:10 2026

@author: goirand-lopez
"""
import os
from os.path import join
import sys
import time

def progressbar(it, prefix="", size=60, out=sys.stdout): # Python3.6+
    print('\n')
    count = len(it)
    start = time.time()
    def show(j):
        x = int(size*j/count)
        remaining = ((time.time() - start) / j) * (count - j)

        mins, sec = divmod(remaining, 60)
        time_str = f"{int(mins):02}:{sec:05.2f}"

        print(f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count} Temps estimé : {time_str}", end='\r', file=out, flush=False)

    for i, item in enumerate(it):
        yield item
        show(i+1)
    print('\n', flush=True, file=out)

# !!! Before using the code you have to Select the right folders!
folderpathname = 'C:/Users/goirand-lopez/Desktop/Neurones-models/'
                
model_folder = join('NEURON 7.8 AMD64', 'DG3D500GCsLFPMeasuresTemp')

cellpath = join(folderpathname,model_folder)
scriptpath = join(folderpathname,join('NEURON 7.8 AMD64', 'DG3D500GCsLFPMeasuresTemp')) # Where your script is
os.chdir(cellpath)
folderpathname = cellpath #'C:/Users/goirand-lopez/Desktop'


PPint = 0
trial = 0

for nspr in [10] :
    
    ParametersSimulation = dict(
    nGClayer = 25,
    PPintervals = PPint ,
    PPnspk = 10 ,
    nprocess = 5,
    simvar = False,
    nspr = nspr,
    lightstim = False,
    GCstim = .01/5, # number of GCs stimulated by each PP,
    IPlist = [0,1,2,3,4],
    randseed = False
    # trial = trial # Variable used to know if we add to load the channel mecanism or if it's already done
    )
    
    with open(os.path.join(cellpath, 'ParametersSimulation'+'.txt'), "w") as f:
        for key,val in ParametersSimulation.items() :
            f.write('{}\t{}\t{}\t'.format(key,val,type(val)))
            f.write('\n')
    
    exec(open(join(scriptpath,'Network500PPs_1IPfor1window_1PP1GC_multiprocess.py')).read())        
    print(
        ''' ---------------------------------------------------
        ********* Waiting to restore cache memory *********
        --------------------------------------------------- ''')
    time.sleep(2*60)