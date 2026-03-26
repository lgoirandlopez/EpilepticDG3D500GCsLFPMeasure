# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 10:52:50 2023
Last Update : Thu Mar 26 2026

@author: goirand-lopez
"""

'''
Here if we have n windows of stimulations, different IPs will be activated within the window
'''    

'''
---------------------------------------------------------------------------------------------------------
*************************************** Librairies import ***************************************
---------------------------------------------------------------------------------------------------------
'''
import multiprocessing
from multiprocessing import Pool, freeze_support, active_children
import LFPy
from LFPy import NetworkCell, Network, NetworkPopulation, Synapse, RecExtElectrode, \
    CurrentDipoleMoment
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
import os
import neuron
from neuron import h
from os.path import join
from mpi4py import MPI
import h5py
import time
import pickle
import datetime
import warnings

from scipy.optimize import minimize
from scipy import integrate as intg
import scipy.sparse as ssparse
import scipy.stats as st
import scipy.signal as ss
import scipy as sc

'''
---------------------------------------------------------------------------------------------------------
*************************************** Function to connect the network cells ***************************************
---------------------------------------------------------------------------------------------------------
'''
def randnorm(binf,bsup,scale=None,size=1) :
    
    loc = (bsup+binf)/2
    if scale is None :
        #' AUTOSCALE'
        scale = (bsup-binf)/6
        
    if size > 1 :
        temp = np.random.normal(loc=loc,scale=scale,size=size)
        temp[temp>bsup]=bsup
        temp[temp<binf]=binf
    else :
        temp = np.random.normal(loc=loc,scale=scale)
        if temp > bsup : 
            temp = bsup
        
        if temp <binf :
            temp = binf
    return temp
        
# *****************************************
#               Connectivity function
# *****************************************
def get_connectivity(network, pretype, posttype, npresyncell, npostsyncell, networkconnparam, 
                     networktype='topographic',randomtype='uniform',ringoption=True,
                     scalePathologicalNetworkRand=0.05) :
    '''
    Function to connect cells of the network according to the parameters of each type of connection.
    
    Parameters
    ----------
    network : Class type
        Distributed populations of cells of type Cell and handling connections between cells in the respective populations.
    pretype : str
        Type of cell of the presynaptic cell.
    posttype : str
        Type of cell of the postsynaptic cell.
    npresyncell : int
        Number of presynaptic cells.
    npostsyncell : int
        Number of possynaptic cells.
    networkconnparam : List of dict
        The list is indexed as followed. 
        Each list correspond to a presynaptic type of cell and within a given list there are as many list as the number of postsynaptic cell type.
        Each list corresponding to a postsynaptic cell type is composed of a dictionnary with every connection parameters.
        The parameters are : 
            * nbconn : the number of connections
            * binf : the inferior born of index cell to pick
            * bsup : the superior born of index cell to pick
    networktype : str, optional
        Type of organisation of the network. The script allows topographic or random.
        The default is 'topographic'.
    randomtype : str, optional
        Type of random used to pick postsynaptic target cell. The default is 'uniform'.
    ringoption : Bool, optional
        Choice to use or not the ring structure for the network. The default is True.
    scalePathologicalNetworkRand : float, optional
        Proportion of random onnection in the topographic recurrent mossy fiber network. The default is 0.05.

    Raises
    ------
    NameError
        Raise an error if an option does not correspond to a waiting one.

    Returns
    -------
    connectivity : Bool matrix
        A boolean matrix of the connections with npresynaptic cell lines and n postsynaptic cell column.

    '''
    nbconn = networkconnparam['nbconn']
    if randomtype.lower() == 'uniform' :
        randfun = np.random.uniform
    elif randomtype.lower() == 'normal' :
        randfun = randnorm
    else :
        raise NameError(randomtype+' is not a right randomtype')
        
    connectivity = np.zeros((npresyncell,npostsyncell))            
    if networktype.lower() == 'topographic' :
      binf = networkconnparam['binf']
      bsup = networkconnparam['bsup']
    elif networktype.lower() == 'random':
        binf = 0
        bsup = npostsyncell
    else :
        raise NameError(networktype+' is not a right networktype')
    
    for npre in range(npresyncell) :
        nconn = 0
        connlist = []
        while nconn < nbconn[npre] :
            if pretype.upper()=='GC' and posttype.upper()=='GC' :
                pick = np.random.rand()
                if pick < scalePathologicalNetworkRand :
                    #'Random connection'
                    shift = int(np.rint(randfun(0,npostsyncell)))  - npre
                else :
                    #'Connection within the topographic window'
                    shift = int(np.rint(randfun(binf,bsup)))
            else :
                shift = int(np.rint(randfun(binf,bsup))) 
            if ('bimodal' in networkconnparam and networkconnparam['bimodal']) :
                shift*=np.random.choice([-1,1])
            if pretype.upper() == posttype.upper() :
                npost = npre + shift
            elif posttype.upper() == 'GC' :
                npost = int(npostsyncell*(npre+.5)//npresyncell) + shift
            elif posttype.upper() == 'MC' :
                if pretype.upper() == 'GC' :
                    npost = int((npostsyncell/3*npre//npresyncell)*3) + shift
                else : 
                    npost = int((npre+1)*npostsyncell//npresyncell) + shift
            elif posttype.upper() == 'BC' :
                npost = int(npostsyncell*npre//npresyncell) + shift
            elif posttype.upper() == 'HC' :
                npost = int(npostsyncell*npre//npresyncell) + shift
            elif posttype.upper() == 'PP' :
                npost = shift #Useful for the case PP for exemple
            if ringoption :
                npost = npost%npostsyncell 
            else : 
                if npost < 0 :
                    npost = 0
                elif npost > 0 :
                    npost = npostsyncell
            if (not ((pretype.upper() == posttype.upper()) and (npost == npre))) and (npost not in connlist) :
                # We avoid autoapse connection and new connection (Il faut ajouter la convergence max)
                #Create the network connection
                connectivity[npre,npost] = 1
                nconn+=1
                connlist.append(npost)
        return connectivity.astype(bool)


def get_stim_Network(network, StimParameters, Stimsynparameters, StimNetworkConnParameters, deadcell = []) :
    '''
    Create the stimulation network and the associated stimulation spike train for each stimulation. 
    The number of different stimulation is set by the NbStimbox.
    Parameters
    ----------
    network : Class type
        Distributed populations of cells of type Cell and handling connections between cells in the respective populations.
    StimParameters : Dict
        General parameters of the stimulation. The different parameters are : 
            * NbStimbox : the number of stimulation channels. 
            * stimstart : the date at which the stimulation begins
            * stimstop : the date at which the stimulation ends
            * RandRegularWindow : the size of the window to randomly shuft the regular spike trains.
            * dt : the time step of the experiment
            * StimTrainType : type of stimulation spike train. Could be regular or random
    Stimsynparameters : List of dict
        List with dict of the synaptic parameters of stimulation for every type of cell to stimulate. The different parameters of each dict are :
            * celltype, the type of cell to connect
            * syntype, the type of synapse for example a double exponential
            * tau1, the rise time (ms)
            * tau2, the decay time (ms)
            * weight, the conductance of the synapse (nS)
            * section, the section of the cell to connect
    StimNetworkConnParameters : List of dict
        List with dict of the network parameters of stimulation for every type of cell to stimulate. The different parameters of each dict are :
            * ConnDiv : Maximal number of presynaptic connections to a target cell from a single stimuli channel
            * nsyn : number of synapste of a connection
            * secreplace : boolean option if a new synapse can replace an older one
            * NetworkType : Two choices - unitary or random. 
                            Unitary means every postsynaptic cell will be connected by a single stimuli channel.
                            Random : Connect randomly each stimulation channel with a target cell according to the set divergence.
    deadcell : List, optional
        List of index of dead cells. The default is [].

    Returns
    -------
    StimConnectivity : List of dict
        A list of dict composed of the name of the target postsynaptic type of cell and the connectivity matrix. 
        The connectivity matrix is a boolean matrix of the connections with npresynaptic cell lines and n postsynaptic cell column.
    SynStimlist : List of list of LFPy synapse class
        List for each stimbox and for each type of target posynaptic cell with a list of synaptic class and its spike times.
    spktrain : Array
        Array of every spike train for each stimulation of the stimbox.

    '''
    SynStimlist = [ [ [] for pp in range(StimParameters['NbStimbox']) ]  for celltype in range(len(Stimsynparameters))]
    StimConnectivity = [dict(celltype=ppsynparam['celltype'],connectivity=[]) for ppsynparam in Stimsynparameters]
    RandRegularWindow = StimParameters['RandRegularWindow']
    StimIntervals = StimParameters['StimIntervals']
    dt = StimParameters['dt']
    stimstart = StimParameters['stimstart']
    stimstop = StimParameters['stimstop']
    # Creation of the stimulus spike trains
    assert StimParameters['StimTrainType'].lower() in ['regular', 'random'], 'Temporal Stimulation type not supported' #
    if StimParameters['StimTrainType'].lower() == 'regular' :
        if RandRegularWindow > 0 and StimIntervals>0:
            # Draw random spikes date within the random regular window if the window exists. 
            # There are n spikes the same for every train and regular window.
            spktrain = np.random.choice(np.arange(0,RandRegularWindow,dt),size=(StimParameters['NbStimbox'],StimParameters['nbspk']))+ np.arange(0,(stimstop-stimstart)/dt,StimIntervals/dt)*dt + stimstart
            spktrain = [s[s<=stimstop] for s in spktrain]
        elif StimIntervals>0 :
            # Draw regular time spaced spikes for every stimulation channels
            assert int((stimstop-stimstart)/StimIntervals)==StimParameters['nbspk'], 'The number of spikes and the regular interval does not correspond'
            spktrain = np.zeros((StimParameters['NbStimbox'],StimParameters['nbspk'])) + np.arange(0,(stimstop-stimstart),StimIntervals)*dt + stimstart
            spktrain = [s[s<=stimstop] for s in spktrain]
        else : 
            # Case of no stimulation window
            spktrain = np.zeros((StimParameters['NbStimbox'],StimParameters['nbspk']))+stimstart
    elif StimParameters['StimTrainType'].lower() == 'random' :
        assert StimIntervals>0, 'Need a mean interval greater than 0'
        # Draw random spikes with a frequency based on the stimulus Interval given
        rnd = np.random.rand(StimParameters['NbStimbox'], int((stimstop-stimstart)/dt))
        lambdat = dt/StimIntervals # scale as in a Poisson process
        spktrain = [np.where(rnd[kk,:]<lambdat)[0]*dt+stimstart for kk in range(StimParameters['NbStimbox'])]
    
    # Creation the stimulus network
    for i, networkconnparam in enumerate(StimNetworkConnParameters) :
        post = Stimsynparameters[i]['celltype']
        npostsyncell = network.populations[post].POP_SIZE
        
        assert networkconnparam['NetworkType'].lower() in ['unitary', 'random'], 'Network type not supported' #
        if networkconnparam['NetworkType'].lower()=='unitary' :
            assert StimParameters['NbStimbox']>=npostsyncell, 'The number of Stimbox as to be less than the number of post synaptic cell, select which to activate after'
            PPDiv = np.arange(npostsyncell) 
        elif networkconnparam['NetworkType'].lower()=='random'  : 
            PPDiv = np.array([np.random.choice(np.arange(0,StimParameters['NbStimbox']),networkconnparam['ConnDiv'],replace=False) for i in range(npostsyncell)])
        connectivity = np.zeros((StimParameters['NbStimbox'],npostsyncell))
        for nc in range(npostsyncell) :
            connectivity[PPDiv[nc],nc] = 1
        if len(deadcell[i]) > 0 :
            connectivity[:,deadcell[i]] = 0
        StimConnectivity[i]['connectivity'] = connectivity
        [idxPP, idxpostsyn] = np.where(connectivity == True)
        for npre,npost in zip(idxPP, idxpostsyn):
            if all(s is None for s in Stimsynparameters[i]['section']) or networkconnparam['nsyn']==0 :
               if all(s is None for s in Stimsynparameters[i]['section']) and not networkconnparam['nsyn']==0 :
                   warnings.warn("You want to connect a synapse to a None section of {} cell".format(Stimsynparameters[i]['celltype']))
               continue
            if len(Stimsynparameters[i]['section'])==networkconnparam['nsyn'] and not networkconnparam['secreplace'] :
                idxbase = network.populations[post].cells[npost].get_idx(Stimsynparameters[i]['section'])
            else :
                idxbase = np.random.choice(network.populations[post].cells[npost].get_idx(Stimsynparameters[i]['section']),
                                       networkconnparam['nsyn'], replace=networkconnparam['secreplace'])        
            for idxpicked in idxbase :
                SynStimlist[i][npre].append(LFPy.Synapse(network.populations[post].cells[npost],
                                              idx=idxpicked,
                                              syntype='Exp2Syn',
                                              weight=Stimsynparameters[i]['weight'],
                                              e=Stimsynparameters[i]['e'],
                                              tau1=Stimsynparameters[i]['tau1'],
                                              tau2=Stimsynparameters[i]['tau2'],                                          
                    ))
                SynStimlist[i][npre][-1].set_spike_times(spktrain[npre])
    
    return StimConnectivity, SynStimlist, spktrain


def HilarLoss(hilarloss,nmcell,nhcell) :
    '''
    Function to kill a fixed rate of cells of the hilus

    Parameters
    ----------
    hilarloss : float
        Rate of loss.
    nmcell : int
        Number of mossy cell.
    nhcell : int
        Number of HIPP.

    Returns
    -------
    DeadMCList : List
        List with the index of killed cells.
    DeadHCList : List
        List with the index of killed cells.

    '''
    DeadMCList = np.random.choice(np.arange(nmcell),int(hilarloss*nmcell),replace=False)
    DeadHCList = np.random.choice(np.arange(nhcell),int(hilarloss*nhcell),replace=False)
    return DeadMCList, DeadHCList

def MCSpontaneousActivity(network, DeadMCList=[]) :
    '''
    Function to simulate spontaneous activity of the mossy cells

    Parameters
    ----------
    network : Class type
        Distributed populations of cells of type Cell and handling connections between cells in the respective populations.
    DeadMCList : List, optional
        List with the index of killed cells. The default is [].

    Returns
    -------
    MCSpontActivity : List 
        List of synapse to activate the MCs with a designed random train spike. 
        Each synapse is followed by a dict with the corresponding information

    '''
    MCSpontActivity = []
    MCSpontActivitySyn = []
    for npost in range(network.populations['MC'].POP_SIZE) :    
        if npost not in DeadMCList :  # don't activate dead cell
            actvscale = 1
        else :
            actvscale = 0
        MCSpontActivitySyn.append(LFPy.Synapse(network.populations['MC'].cells[npost],
                           idx = 0, #soma
                            syntype='Exp2Syn',
                            weight=actvscale*10e-3, #1.275e-3,
                            e=0,
                            tau1=1.5,
                            tau2=5.5))
        # syn.set_spike_times(np.array([200]))
        freqMC = 0
        while freqMC<2 or freqMC>4 :
            freqMC = np.random.normal(3,scale=.5) # random freq between 2-4 HZ (si au dessus ou en dessous on relance)
        
        noise = .75
        nspkMC = int(np.random.normal(freqMC*network.tstop//1000))
        MCSpontActivitySyn[-1].set_spike_times_w_netstim(noise=noise, start=0, number=nspkMC, interval=1000/freqMC)
        MCSpontActivity.append(dict(freqMC=freqMC,nspkMC=nspkMC,noise=noise))
    return MCSpontActivity

def expfunsyn(gsyn,v,t,tsp,delay,tauon,tauoff,erev):
    '''
    Function to simulate a synaptic current. The synapptic current is the two state kinetic scheme synapse described by rise time tau1,
    and decay time constant tau2. The normalized peak condunductance is 1.
    Decay time MUST be greater than rise time.

    Parameters
    ----------
    gsyn : float
        Conductance of the synapse (μS).
    v : float
        The voltage at the moment t (mV).
    t : float
        The timing moment (ms).
    tsp : float
        The time of presynaptic cell spike (ms).
    delay : float
        The timing delay due to conduction between the presynaptic spike and the posynaptic current (ms).
    tauon : float
        Rise time (ms).
    tauoff : float
        Decay time (ms)
    erev : float
        The reversal potential of the synapse (mV).

    Returns
    -------
    Current : float
        Postsynaptic current at the moment t (nA).

    '''
    tp = (tauon*tauoff)/(tauoff - tauon) * np.log(tauoff/tauon)
    fnorm = 1/(np.exp(-tp/tauoff) - np.exp(-tp/tauon))
    return gsyn*fnorm*(np.exp(-(t-tsp-delay)/tauoff)-np.exp(-(t-tsp-delay)/tauon))*(v-erev)*(t>tsp)

'''
---------------------------------------------------------------------------------------------------------
***************************************** Experiment Parameters *****************************************
---------------------------------------------------------------------------------------------------------
'''

    

def custom_callback(result):
    '''
    Function to give the result of the simulation

    Parameters
    ----------
    result : str
        The message displays when a simulation with the parallel process is over.

    Returns
    -------
    None.

    '''
    print(f'Got result: {result}')

def custom_error_callback(error):
 	print(f'Got error: {error}')
        

#!!! Main code
def Network500PPs_1IPfor1window_1PP1GC_PPintervalExploration(IP, nspr, PPintervals, KARfactor, PPnspk,
                                                             resultsfolderpathname, cellpath,
                                                             GCstim=1, nGClayer=25, simvar=True, 
                                                             lightstim=False, randseed=False, saveopt=True) :
    '''
    Simulation function of the reduced dentate gyrus as introduced in the article.

    Parameters
    ----------
    IP : int
        Identification number of the simulation.
    nspr : int
        The sprouting of the pathological network. A sprouting of n means that each GC will connect n other GCs.
    PPintervals : int
        Intervals between spike og the PP stimulations (ms).
    KARfactor : float
        Rate of kainate receptors, between 0 and 1.
    PPnspk : int
        Number of spikes in the PP stimulation.
    resultsfolderpathname : str
        Path of the folder where the results will be stored.
    cellpath : str
        Path of the folder where the cells models are.
    GCstim : Int, optional
        Number of GCs stimulate by each PP stimulation (PP stim box). The default is 1.
    nGClayer : int
        Number of Gcs by layer of the simulated dentate gyrus. The default is 25.
    simvar : Bool, optional
        Boolean option to simulate the activity or not. The default is True.
    lightstim : Bool, optional
        Boolean option to add a light stimulation during the simulation. Useful to test the reaction of pseudostationnarity activity. 
        The default is False.
    randseed : Bool, optional
        Boolean option to fixe a new random seed or not.
    saveopt : Bool, optional
        Boolean option to save or not the results.

    Returns
    -------
    None.

    '''
    
    COMM = MPI.COMM_WORLD
    RANK = int(randseed)*COMM.Get_rank() + int(IP+(nspr+1)*(PPnspk)) 
    # If not randseed : cancel Rank to fix the seed
    # if randeed : 
    # avoid same sequence of random numbers from numpy and neuron on each RANK,
    # e.g., in order to draw unique cell and synapse locations and random synapse
    # activation times
    GLOBALSEED = 1234*(1-int(randseed)) + int(randseed)*(int(time.time())-1000*int(time.time()//1000))
    # If not randseed : Fixed Seed at 1234, if randseed 
    np.random.seed(GLOBALSEED + RANK) # select the random seed for every random element
    
    
    '''
----------------------------------------------------------------------------------------------
***************************************** Set the global parameters *****************************************
----------------------------------------------------------------------------------------------
'''
    ngcell = 500
    nmcell = 15
    nbcell = 6
    nhcell = 6
    
    dt = 2**-3
    tstart = - 150 #500 
    
    # Base for the filename and the specific result path and folder
    Filenametemp = 'ResultatsNetworkExperiences'
    if Filenametemp not in os.listdir(resultsfolderpathname) :
            OUTPUTPATH=join(resultsfolderpathname,Filenametemp)
            os.mkdir(OUTPUTPATH)
    BASEPATH4RESULT = join(resultsfolderpathname,Filenametemp)
    
    Filenametemp = '{:n}Layers{:n}GCbyLayer'.format(ngcell/nGClayer,nGClayer)
    if Filenametemp not in os.listdir(BASEPATH4RESULT) :
            OUTPUTPATH=join(BASEPATH4RESULT,Filenametemp)
            os.mkdir(OUTPUTPATH)
    BASEPATH4RESULT = join(BASEPATH4RESULT,Filenametemp)    
    
    '''
----------------------------------------------------------------------------------------------
***************************************** Set the format of image *****************************************
----------------------------------------------------------------------------------------------
'''
    dpi=300
    imgfrm = '.svg' #'.eps'

    
    '''
----------------------------------------------------------------------------------------------
***************************************** Set the network options *****************************************
----------------------------------------------------------------------------------------------
'''
    InhibitionFactor = 1 # scaling of synaptic weight BC->X and HC->X
    IFLFactor = 1 #3 #InhibitoryFeedbackLoopFactor  # scaling of synaptic weight GC->BC and BC->GC
    
    scale_GC2GC = 4     # scaling of synaptic weight GC->GC 
    scale_MC2GC = 1     # scaling of synaptic weight MC->GC 
    scale_BC2GC = IFLFactor    # scaling of synaptic weight BC->GC 
    scale_HC2GC = IFLFactor      # scaling of synaptic weight HC->GC (beta_HIPP in Myers and Scharfman, 2009)
    
    SpontMCFactor = 0
    
    
    KARopt = KARfactor>0 #Presence or absence of KARs for tilte purpose
    scaleKAR = .33 # Scale of synaptic weight of KAR synapse compare to AMPAR one
    NaPopt = True # Activate NaP or not in the GCs
    
    hilarloss = 1*.5 # Rate of hilar loss
    networktype = 'Topographic' # Type of network can be 'Topographic' or 'random'
    ringoption = True # Option or ring structure for network
    scalePathologicalNetworkRand = 0.05
    randomtype='Uniform' # Type of random used for the selection of target cells within a range
    
    NetworkParameters2 = dict(InhibitionFactor=InhibitionFactor,InhibitoryFeedbackLoopFactor=IFLFactor,
                              SpontaneousMCFactor=SpontMCFactor,sprouting=nspr,NaPOption=NaPopt,RingOption=ringoption,
                              NetworkType=networktype, RandomType=randomtype,
                              scalePathologicalNetworkRand=scalePathologicalNetworkRand,
                              scale_BC2GC=scale_BC2GC,scale_HC2GC=scale_HC2GC,scale_MC2GC=scale_MC2GC)
    
    '''
----------------------------------------------------------------------------------------------
***************************************** Set the PPStim options *****************************************
----------------------------------------------------------------------------------------------
'''
    npp = 500 # Number of GCs PP
    stimstart = 10
    stimstop =  PPnspk*PPintervals+stimstart #tstop #
    tstop = stimstop + 300 #stimstop + 200
    randstimwind = 0 # Size of the window to randomly shift a regular stimulus spike train
    
    # scale of stiumulation a scale of 0 means no stimulation
    scale_PP2GC = 1
    scale_PP2MC = 0
    scale_PP2BC = 1
    scale_PP2HC = 0
    
    PPTraintype='Regular'
    
    '''
----------------------------------------------------------------------------------------------
******************************** Results Name Definition and pathname definition ********************************
----------------------------------------------------------------------------------------------
'''
    if ringoption :
        ringname = 'Ring'
    else : 
        ringname = 'NoRing'
        
    if (InhibitionFactor>0) :
        inhibitionName = "Inhibited"
    else :
    	inhibitionName = "Disinhibited"
    
    ResultFile = 'Scale[GC-GC, MC-GC, BC-GC, HC-GC]={}'.format([scale_GC2GC,scale_MC2GC,
                                                              scale_BC2GC,scale_HC2GC])
    if lightstim :
        ResultFile2 = 'SparseInputsRegularIntervalAndLightStim' 
    else :
        ResultFile2 = 'SparseInputsRegularInterval'
    
    
    if hilarloss > 0 :
        ResultFileSclerosis = 'Sclerosis'
    else :
        ResultFileSclerosis =  'No sclerosis'
        
    
    if nspr == 0 :
        if 'HealthyDG' not in os.listdir(BASEPATH4RESULT) :
            OUTPUTPATH=join(BASEPATH4RESULT,'HealthyDG')
            os.mkdir(OUTPUTPATH)
        else :
            OUTPUTPATH=join(BASEPATH4RESULT,'HealthyDG')
        
        if ResultFile not in os.listdir(OUTPUTPATH) :
            OUTPUTPATH=join(OUTPUTPATH,ResultFile)
            os.mkdir(OUTPUTPATH)
        else :
            OUTPUTPATH=join(OUTPUTPATH,ResultFile)
        
        if ResultFile2 not in os.listdir(OUTPUTPATH) :
            OUTPUTPATH=join(OUTPUTPATH,ResultFile2)
            os.mkdir(OUTPUTPATH)
            os.mkdir(join(OUTPUTPATH,'Data'))
        else :
            OUTPUTPATH=join(OUTPUTPATH,ResultFile2)
        
        
        if NaPopt :
            basename = "network_Healthy_NaP"  
        else : 
            basename = "network_Healthy_noNaP"  
    
                                                     
    elif KARopt :
        if 'EpilepticDGNaPKaR'+'_scaleKA{}%'.format(int(100*scaleKAR)) not in os.listdir(BASEPATH4RESULT) :
            OUTPUTPATH=join(BASEPATH4RESULT,'EpilepticDGNaPKaR'+'_scaleKA{}%'.format(int(100*scaleKAR)))
            os.mkdir(OUTPUTPATH)
        else :
            OUTPUTPATH=join(BASEPATH4RESULT,'EpilepticDGNaPKaR'+'_scaleKA{}%'.format(int(100*scaleKAR)))
            
        if 'EpilepticDGNaPKaR'+'_scaleKA{}%'.format(int(100*scaleKAR)) not in os.listdir(BASEPATH4RESULT) : 
            os.mkdir(OUTPUTPATH)
        if ResultFile not in os.listdir(OUTPUTPATH) :
            OUTPUTPATH=join(OUTPUTPATH,ResultFile)
            os.mkdir(OUTPUTPATH)
        else :
            OUTPUTPATH=join(OUTPUTPATH,ResultFile)
        
        if ResultFileSclerosis not in os.listdir(OUTPUTPATH) :
            OUTPUTPATH=join(OUTPUTPATH,ResultFileSclerosis)
            os.mkdir(OUTPUTPATH)
        else :
            OUTPUTPATH=join(OUTPUTPATH,ResultFileSclerosis)
        
        if ResultFile2 not in os.listdir(OUTPUTPATH) :
            OUTPUTPATH=join(OUTPUTPATH,ResultFile2)
            os.mkdir(OUTPUTPATH)
            os.mkdir(join(OUTPUTPATH,'Data'))
        else :
            OUTPUTPATH=join(OUTPUTPATH,ResultFile2)
            
    
        if NaPopt :
            basename = "network_KAR_NaP-{}%scaleKA".format(100*KARfactor)
        else :
            basename = "network_KAR_noNaP-{}%scaleKA".format(100*KARfactor)
            
    else :
        if 'EpilepticDGNaPAMPAROnly' not in os.listdir(BASEPATH4RESULT) :
            OUTPUTPATH=join(BASEPATH4RESULT,'EpilepticDGNaPAMPAROnly')
            os.mkdir(OUTPUTPATH)
        else :
            OUTPUTPATH=join(BASEPATH4RESULT,'EpilepticDGNaPAMPAROnly')
        
        if ResultFile not in os.listdir(OUTPUTPATH) :
            OUTPUTPATH=join(OUTPUTPATH,ResultFile)
            os.mkdir(OUTPUTPATH)
        else :
            OUTPUTPATH=join(OUTPUTPATH,ResultFile)
        
        if ResultFileSclerosis not in os.listdir(OUTPUTPATH) :
            OUTPUTPATH=join(OUTPUTPATH,ResultFileSclerosis)
            os.mkdir(OUTPUTPATH)
        else :
            OUTPUTPATH=join(OUTPUTPATH,ResultFileSclerosis)
        
        if ResultFile2 not in os.listdir(OUTPUTPATH) :
            OUTPUTPATH=join(OUTPUTPATH,ResultFile2)
            os.mkdir(OUTPUTPATH)
            os.mkdir(join(OUTPUTPATH,'Data'))
        else :
            OUTPUTPATH=join(OUTPUTPATH,ResultFile2)
    
        if NaPopt :
            basename = "network_AMPAR_NaP"
        else :
            basename = "network_AMPAR_noNaP"
    
    if lightstim :
        idname = "%s_%s-%d%s-%d%s-%d%s-%d%s-%d%s-%d%s-%d%s-%s" %(ringname+networktype+randomtype+"Random{}%".format(int(100*scalePathologicalNetworkRand))+"PPTrainRand"+inhibitionName,
                                                    basename,nspr,"%sprounting",int(hilarloss*100),'%HilarCellLoss',
                                                    GCstim,"GCsStim",
                                                    PPintervals,"PPinterval",PPnspk,"PPnspk",
                                                    randstimwind,"RandStimWindow",int(scale_PP2GC*100),"%wPP2GC",'LightStim')
    else : 
        idname = "%s_%s-%d%s-%d%s-%d%s-%d%s-%d%s-%d%s-%d%s-%s" %(ringname+networktype+randomtype+"Random{}%".format(int(100*scalePathologicalNetworkRand))+"PPTrainRand"+inhibitionName,
                                                basename,nspr,"%sprounting",int(hilarloss*100),'%HilarCellLoss',
                                                GCstim,"GCsStim",
                                                PPintervals,"PPinterval",PPnspk,"PPnspk",
                                                randstimwind,"RandStimWindow",int(scale_PP2GC*100),"%wPP2GC",'Base')
       
    idnameIP = '-IP{}-1PP1GC'.format(IP)+idname
    OUTPUTPATH2=join(OUTPUTPATH,'Data')
    
                        
    '''
    ---------------------------------------------------------------------------------------------------------
    *********************************** Set cells parameters ************************************
    ---------------------------------------------------------------------------------------------------------
    '''
    
    networkParameters = dict(
        dt =dt,
        tstart = tstart,
        tstop = tstop, #1500
        v_init = -70.,
        celsius = 6.3,
        OUTPUTPATH = OUTPUTPATH2,
            )
    
    GC_parameters = {
                    'morphology': join('GC0_ep.hoc'),
                    'templatefile': join('GCtemplate.hoc'),
                    'templatename': 'GranuleCell',
                    'templateargs': None,
                    'v_init': -80,
                    'passive': False,
                    'nsegs_method': None,
                    'delete_sections': False,
                    'pt3d': True,
                    'tstart': tstart,  # [ms] Simulation start time
                    'dt': dt,  # [ms] Should be a power of 2
                    'tstop': tstop,  # [ms] Simulation end time
            }
    if  nspr > 0 : 
        if NaPopt :
            GC_parameters['morphology'] = join('GC0_ep.hoc')
        else : 
            GC_parameters['morphology'] = join('GC0_ep_noNaP.hoc')
    else : 
        if NaPopt :
            GC_parameters['morphology'] = join('GC0.hoc')
        else : 
            GC_parameters['morphology'] = join('GC0_noNaP.hoc')
        
    
    MC_parameters = {
                    'morphology': join('MC.hoc'),
                    'templatefile': join('MCtemplate.hoc'),
                    'templatename': 'MossyCell',
                    'templateargs': None,
                    'v_init': -60,
                    'passive': False,
                    'nsegs_method': None,
                    'delete_sections': False,
                    'pt3d': True,
                    'dt': dt,  # [ms] Should be a power of 2
                    'tstart': tstart,  # [ms] Simulation start time
                    'tstop': tstop,  # [ms] Simulation end time
            }
    
    BC_parameters = {
                    'morphology': join('BC.hoc'),
                    'templatefile': join('BCtemplate.hoc'),
                    'templatename': 'BasketCell',
                    'templateargs': None,
                    'v_init': -60,
                    'passive': False,
                    'nsegs_method': None,
                    'delete_sections': False,
                    'pt3d': True,
                    'dt': dt,  # [ms] Should be a power of 2
                    'tstart': tstart,  # [ms] Simulation start time
                    'tstop': tstop,  # [ms] Simulation end time
            }
    
    HC_parameters = {
                    'morphology': join('HC.hoc'),
                    'templatefile': join('HCtemplate.hoc'),
                    'templatename': 'HIPPCell',
                    'templateargs': None,
                    'v_init': -70,
                    'passive': False,
                    'nsegs_method': None,
                    'delete_sections': False,
                    'pt3d': True,
                    'dt': dt,  # [ms] Should be a power of 2
                    'tstart': tstart,  # [ms] Simulation start time
                    'tstop': tstop,  # [ms] Simulation end time
            }
    
    
    GCpopulationParameters = dict(
        Cell=NetworkCell,
        cell_args=GC_parameters,
        pop_args=dict(
            radius=100.,
            loc=0.,
            scale=20.),
        rotation_args=dict(x=0., y=0.),
    )
    
    MCpopulationParameters = dict(
        Cell=NetworkCell,
        cell_args=MC_parameters,
        pop_args=dict(
            radius=100.,
            loc=0.,
            scale=20.),
        rotation_args=dict(x=0., y=0.),
    )
    
    BCpopulationParameters = dict(
        Cell=NetworkCell,
        cell_args=BC_parameters,
        pop_args=dict(
            radius=100.,
            loc=0.,
            scale=20.),
        rotation_args=dict(x=0., y=0.),
    )
    
    HCpopulationParameters = dict(
        Cell=NetworkCell,
        cell_args=HC_parameters,
        pop_args=dict(
            radius=100.,
            loc=0.,
            scale=20.),
        rotation_args=dict(x=0., y=0.),
    )
    
    
    population_names = ['GC', 'MC', 'BC', 'HC']
    population_sizes = [ngcell, nmcell, nbcell, nhcell]
    CelltypeParameters = [GC_parameters, MC_parameters, BC_parameters, HC_parameters]
    populationParameters = [GCpopulationParameters, MCpopulationParameters, BCpopulationParameters, HCpopulationParameters]
    
    '''
    ---------------------------------------------------------------------------------------------------------
    ******************************** Stimulation Parameters (PP Stimulation) ********************************
    ---------------------------------------------------------------------------------------------------------
    '''
    # ******************************************************
    #          Connectivity and Synapse Parameters
    # ******************************************************
    
    PPsynParameters = [dict(celltype='GC', syntype=neuron.h.Exp2Syn, tau1=1.5, tau2=5.5, e=0, weight = scale_PP2GC*20e-3,
                            section=['gcdend1[3]', 'gcdend2[3]']),
                       dict(celltype='MC',syntype=neuron.h.Exp2Syn, tau1=1.5, tau2=5.5, e=0, weight = scale_PP2MC*5e-3,
                            section=['mcdend1[3]', 'mcdend2[3]', 'mcdend3[3]', 'mcdend4[3]']),
                       dict(celltype='BC',syntype=neuron.h.Exp2Syn, tau1=2, tau2=6.3, e=0, weight = scale_PP2BC*10e-3,
                            section=['bcdend1[3]', 'bcdend2[3]']),
                       dict(celltype='HC',syntype=neuron.h.Exp2Syn, tau1=1.5, tau2=5.5, e=0, weight = scale_PP2HC*5e-3,
                            section=[None])]
    
    PPStimParameters = dict(NbStimbox=npp, StimIntervals=PPintervals, nbspk=PPnspk, 
                          stimstart=stimstart, RandRegularWindow=randstimwind, stimstop=stimstop, dt=dt,
                          StimTrainType=PPTraintype)
    
    
    PPStimNetworkConnParameters = [dict(ConnDiv=1, nsyn=2, secreplace=False, NetworkType='unitary'),
                              dict(ConnDiv=0, nsyn=1, secreplace=False, NetworkType='random'),
                              dict(ConnDiv=int(npp/6), nsyn=1, secreplace=False, NetworkType='random'),
                              dict(ConnDiv= 0, nsyn=1, secreplace=False, NetworkType='random')]
    
    '''
    ---------------------------------------------------------------------------------------------------------
    ******************************** Light excitatory Stimulation Parameters (Stimulation to inspect chaotic dynamics) ********************************
    ---------------------------------------------------------------------------------------------------------
    '''
    # ******************************************************
    #             Connectivity and Synapse Parameters
    # ******************************************************
    
    LightStimsynParameters = [dict(celltype='GC', syntype=neuron.h.Exp2Syn, tau1=1.5, tau2=5.5, e=0, weight = scale_PP2GC*20e-3,
                            section=['gcdend1[3]', 'gcdend2[3]'])]
    
    LightStimParameters = dict(NbStimbox=npp, StimIntervals=0, nbspk=1, 
                          stimstart=stimstop/2, RandRegularWindow=randstimwind, stimstop=stimstop/2+dt, dt=dt,
                          StimTrainType=PPTraintype)
        
    LightStimNetworkConnParameters = [dict(ConnDiv=1, nsyn=2, secreplace=False, NetworkType='unitary')]
    
     
    
    
    '''
    ---------------------------------------------------------------------------------------------------------
    *********************************** Set Synapstic parameters ***********************************
    ---------------------------------------------------------------------------------------------------------
    '''
    
    # Spike Detection Parameters : allow to choose for each neuronal type the spiking threshold 
    
    #threshold : (mV), delay : (ms)
    SpikeDetectionParam = [dict(target=None,threshold=10.0, weight=0.0, delay=0.0),
                           dict(target=None,threshold=10.0, weight=0.0, delay=0.0),
                           dict(target=None,threshold=-10.0, weight=0.0, delay=0.0),
                           dict(target=None,threshold=10.0, weight=0.0, delay=0.0)]
    
    synapseModel = neuron.h.Exp2Syn
    
    #tau : (ms), e : (mV)
    synapseParameters = [[dict(tau1=1.2, tau2=4, e=0),
                          dict(tau1=.5, tau2=6.2, e=0),
                          dict(tau1=.3, tau2=.5, e=0),
                          dict(tau1=.3, tau2=.6, e=0)],
                         [dict(tau1=1.5, tau2=5.5, e=0),
                          dict(tau1=.45, tau2=2.2, e=0),
                          dict(tau1=.9, tau2=3.6, e=0),
                          dict(tau1=.9, tau2=3.6, e=0)],
                         [dict(tau1=.26, tau2=5.5, e=-70),
                          dict(tau1=.3, tau2=3.3, e=-70),
                          dict(tau1=.16, tau2=1.8, e=-70),
                          dict(tau1=None, tau2=None, e=None)],
                         [dict(tau1=.5, tau2=6, e=-70),
                          dict(tau1=.5, tau2=6, e=-70),
                          dict(tau1=.4, tau2=5.8, e=-70),
                          dict(tau1=None, tau2=None, e=None)]
                         ]
    GCGCsynParameters = [dict(tau1=1.2, tau2=4.25,e=0), dict(tau1=3, tau2=45,e=0)] # [AMPA, KA]
    if scale_GC2GC <1 :
        GCGCweightArguments = [dict(loc=scale_GC2GC*1.75e-3, scale=0), dict(loc=scaleKAR*scale_GC2GC*1.75e-3, scale=0)] # [AMPA, KA]
    else :
        GCGCweightArguments = [dict(loc=scale_GC2GC*1.75e-3, scale=0), dict(loc=scale_GC2GC*scaleKAR*1.75e-3, scale=0)] # [AMPA, KA]
    
    
    # Weight in nS
    weightFunction = np.random.normal
    weightArguments = [[dict(loc=scale_GC2GC*1.75e-3, scale=0),
                        dict(loc=.2e-3, scale=0),
                        dict(loc=4.7e-3, scale=0),
                        dict(loc=.5e-3, scale=0)],
                       [dict(loc=scale_MC2GC*.3e-3, scale=0),
                        dict(loc=.5e-3, scale=0),
                        dict(loc=.3e-3, scale=0),
                        dict(loc=.2e-3, scale=0)],
                       [dict(loc=InhibitionFactor*scale_BC2GC*1.6e-3, scale=0),
                        dict(loc=InhibitionFactor*1.5e-3, scale=0),
                        dict(loc=InhibitionFactor*7.6e-3, scale=0),
                        dict(loc=None, scale=0)],
                       [dict(loc=InhibitionFactor*scale_HC2GC*.5e-3, scale=0),
                        dict(loc=InhibitionFactor*1.5e-3, scale=0),
                        dict(loc=InhibitionFactor*.5e-3, scale=0),
                        dict(loc=None, scale=0)]
                       ]
    minweight = 0
    
    # delay : (ms)
    delayFunction = sc.stats.truncnorm # np.random.normal
    delayArguments = [[dict(a=0, b=np.inf, loc=.8, scale=0),
                       dict(a=0, b=np.inf, loc=1.5, scale=0),
                      dict(a=0, b=np.inf, loc=.8, scale=0),
                       dict(a=0, b=np.inf, loc=1.5, scale=0)],
                      [dict(a=0, b=np.inf, loc=3, scale=0),
                       dict(a=0, b=np.inf, loc=2, scale=0),
                      dict(a=0, b=np.inf, loc=3, scale=0),
                       dict(a=0, b=np.inf, loc=3, scale=0)],
                      [dict(a=0, b=np.inf, loc=.85, scale=0),
                       dict(a=0, b=np.inf, loc=1.5, scale=0),
                      dict(a=0, b=np.inf, loc=.8, scale=0),
                       dict(a=0, b=np.inf, loc=None, scale=0)],
                      [dict(a=0, b=np.inf, loc=1.6, scale=0),
                       dict(a=0, b=np.inf, loc=1, scale=0),
                      dict(a=0, b=np.inf, loc=1.6, scale=0),
                       dict(a=0, b=np.inf, loc=None, scale=0)]
                      ]
    
    mindelay = None # if 0, it will be deprecated

    
    networkConnParameters = [[dict(nbconn=ngcell*[nspr], binf=-50, bsup=50),
                              dict(nbconn=ngcell*[1], binf=-2, bsup=2),
                              dict(nbconn=ngcell*[1], binf=-1, bsup=1),
                              dict(nbconn=ngcell*[3], binf=-2, bsup=2)],
                             [dict(nbconn=nmcell*[200], binf=175, bsup=25, bimodal=True),
                              dict(nbconn=nmcell*[3], binf=-3, bsup=3),
                              dict(nbconn=nmcell*[1], binf=-3, bsup=3),
                              dict(nbconn=nmcell*[2], binf=-2, bsup=2)],
                             [dict(nbconn=nbcell*[100], binf=-70, bsup=70),
                              dict(nbconn=nbcell*[3], binf=-3, bsup=3),
                              dict(nbconn=nbcell*[2], binf=-1, bsup=1),
                              dict(nbconn=nbcell*[0], binf=None, bsup=None)],
                             [dict(nbconn=nhcell*[160], binf=-130, bsup=130),
                              dict(nbconn=nhcell*[4], binf=-2, bsup=2),
                              dict(nbconn=nhcell*[4], binf=-2, bsup=2),
                              dict(nbconn=nhcell*[0], binf=None, bsup=None)]]
    
    synapsePositionArguments = [[dict(section=['gcdend1[1]', 'gcdend2[1]']),
                                 dict(section=['mcdend1[1]', 'mcdend2[1]', 'mcdend3[1]', 'mcdend4[1]']),
                                 dict(section=['bcdend1[0]', 'bcdend2[0]', 'bcdend3[0]', 'bcdend4[0]']),
                                 dict(section=['hcdend1[0]', 'hcdend2[0]', 'hcdend3[0]', 'hcdend4[0]'])],
                                [dict(section=['gcdend1[1]', 'gcdend2[1]']),
                                 dict(section=['mcdend1[0]', 'mcdend2[0]', 'mcdend3[0]', 'mcdend4[0]']),
                                 dict(section=['bcdend1[1]', 'bcdend2[1]']),
                                 dict(section=['hcdend1[1]', 'hcdend2[1]', 'hcdend3[1]', 'hcdend4[1]'])],
                                [dict(section=['soma[0]']),
                                 dict(section=['soma[0]']),
                                 dict(section=['bcdend1[1]', 'bcdend2[1]']),
                                 dict(section=[None])],
                                [dict(section=['gcdend1[3]', 'gcdend2[3]']),
                                 dict(section=['mcdend1[2]', 'mcdend2[2]', 'mcdend3[2]', 'mcdend4[2]']),
                                 dict(section=['bcdend1[3]', 'bcdend2[3]']),
                                 dict(section=[None])]
                                ]
    
      
    
    
        
    print('''
    ---------------------------------------------------------------------------------------------------------
    ******************************************** Initiate Network ********************************************
    ---------------------------------------------------------------------------------------------------------
    ''')
    # ********************************************
    #            Network Object Creations
    # ********************************************
    network = Network(**networkParameters)
    
    
    
    for name, size, popParam in zip(population_names, population_sizes, populationParameters):
        print('------- '+name+' -------')
        network.create_population(name=name, POP_SIZE=size, CELLPATH=cellpath,idname=idnameIP,
                                  **popParam)        
    # *****************************************
    #          Hilar Cell Loss
    # *****************************************
    
    if hilarloss>0 :
        DeadMCList, DeadHCList = HilarLoss(hilarloss,nmcell,nhcell)
        DeadCell = [[],DeadMCList,[],DeadHCList]
    else :
        DeadCell = [[] for i in range(len(population_names))]
    
    print('''
    ---------------------------------------------------------------------------------------------------------
    ******************************************** Set the position of the cells ********************************************
    ---------------------------------------------------------------------------------------------------------
    ''')
    # **************************************************************************
    #               FUNCTIONS to calculate GC rotation for ecah cells
    # **************************************************************************               
    def lengthfun(t) :
        return np.sqrt((r1*np.sin(t))**2+(r2*np.cos(t))**2) # We only keep the elliptic function
    # To obtain an hyperoblic function the cosisnus and sinus has to be replaced by hyperbolic cosinus and sinus
    
    def tmaxresearch(tmax,t0,lengthmax) :
        return abs(lengthmax-intg.quadrature(lengthfun, t0, tmax)[0])
    
    def normvec(r1,r2,t) :
        return np.sqrt((np.cos(t)/r1)**2+(np.sin(t)/r2)**2)
    
    def prodscalaire(r1,r2,t0_,t) :
        return (np.cos(t0_)*np.cos(t)/r1**2 + np.sin(t0_)*np.sin(t)/r2**2)/(normvec(r1,r2,t0_)*normvec(r1,r2,t))
    
    # *****************************************
    #              Septotemporal Axis Parameters
    # *****************************************
    LayerPop = [nGClayer,1,1,1] # Number of cells by layer where the cells are present
    septotemporalextent = 6000 # um Septotemporal extent of the GCL
    
    
    # *****************************************
    #              Elliptic GCL Model Parameters
    # *****************************************
    
    t0 = 0.48*np.pi # Minimal Angle 
    t1 = 1.49*np.pi # Maximal Angle
    tlist = np.linspace(t0,t1,num = 100)
    
    r1 = 826.5 # Semi-Major axis
    r2 = 202.7 # Semi-Minor axis
    
    #Create the GCLline and its borders 
    widthGCL = 80 #um
    borderGCLHilus = []
    borderMLGCL = []
    
    for t in tlist :
        x0 = r1*np.cos(t)
        y0 = r2*np.sin(t)
        # ------ projection on orthonal vector of the points to obtain both side of the borders for the point ------
        w0 = widthGCL/np.sqrt(1+(r1**2/r2**2*y0/x0)**2)/2 # projected width on the x-axis
        # w0/2 because the distance between the GCl and the borders is the width/2
        xHilus = -1*np.sign(np.cos(t))*w0+x0 # 
        yHilus = y0*(r1**2/r2**2*(xHilus-x0)/x0+1)
        borderGCLHilus.append([xHilus,yHilus])
        xML = np.sign(np.cos(t))*w0+x0
        yML = y0*(r1**2/r2**2*(xML-x0)/x0+1)
        borderMLGCL.append([xML,yML])
    
    borderGCLHilus = np.array(borderGCLHilus)
    borderMLGCL = np.array(borderMLGCL)
    
    
    # Length of Simulated GCL
    lengthGCL = intg.quadrature(lengthfun, t0, t1)[0] #Lengthfun with the same mtric than the figure
    print('Length of GCL = {:.2f} um'.format(lengthGCL))
    
    GCinterv = int(10*(lengthGCL/nGClayer//10))
    
    # Angle Calculation for the 5 GCS in GCL, the 5 soma are separated with 300um, the first and last GC are from 150um to end of GCL
    tGClist = []
    tmax = minimize(tmaxresearch,t0+0.5, args=(t0,GCinterv//2),tol=1e-6)
    print('tmax = {}'.format(tmax.x[0]))
    tGClist.append(tmax.x[0])
    for ind in range(nGClayer-1) :
        tmax = minimize(tmaxresearch,tmax.x[0]+np.round((t1-t0)/(nGClayer-1),decimals=2), args=(tmax.x[0],GCinterv),tol=1e-6)
        print('tmax = {}'.format(tmax.x[0]))
        tGClist.append(tmax.x[0])
    
    [xsomamatGC,ysomamatGC] = [r1*np.cos(np.array(tGClist)),r2*np.sin(np.array(tGClist))]
    # Calculation of angle between first GC normal to gCL and others
    t0_ = tGClist[0] #Start Angle (first GC)
    tGClistRotation = [2*np.pi*((t-t0_)>np.pi) + (1-2*((t-t0_)>np.pi)) * np.arccos(np.round(prodscalaire(r1,r2,t0_,t),4)) for t in tGClist]
    
    
    secname = 'Soma'
    geompath = join(cellpath,'CoordBySectionCell')
    
    rotname = os.path.join(OUTPUTPATH2, 'cell_positions_and_rotations'+idnameIP+'.h5')
    posrotfile = h5py.File(rotname,'r')
    
    for cellname, ncell, nbcelllayer in zip(population_names,population_sizes,LayerPop) :   
        print("****************** Geometry of "+ cellname +" ******************")
        zwidthLayer = septotemporalextent/int(ncell/nbcelllayer)
        zLayers = zwidthLayer/2 + zwidthLayer*np.arange(0,int(ncell/nbcelllayer))
        print ("----- " + cellname +" ----- Number of layers : {} - Between each layer you have {}um on the z-axis".format(int(ncell/nbcelllayer),zwidthLayer))
        
        # Load the localisation of the cell
        if cellname!='GC' :
            fd = open(join(geompath,cellname+'x'+secname+'.txt'),mode='r')
            xsomamat = np.loadtxt(fd)
            fd.close()      
            
            fd = open(join(geompath,cellname+'y'+secname+'.txt'),mode='r')
            ysomamat = np.loadtxt(fd)
            fd.close()   
        else : 
            xsomamat = xsomamatGC
            ysomamat = ysomamatGC 
        
        
        for indlayer in range(int(ncell/nbcelllayer)) : 
            for indcelllayer in range(nbcelllayer) :
                indcell = indlayer*nbcelllayer+indcelllayer
                if nbcelllayer>1 :
                    xsoma = [xsomamat[indcelllayer],0]
                    ysoma = [ysomamat[indcelllayer],0]
                else :
                    xsoma = xsomamat
                    ysoma = ysomamat
                network.populations[cellname].cells[indcell].set_pos(x=xsoma[0],y=ysoma[0],z=zLayers[indlayer])
                network.populations[cellname].cells[indcell].set_rotation(x=0,y=0,
                                                                          z=tGClistRotation[indcelllayer]-posrotfile[cellname][indcell][-1])
    posrotfile.close()
        
    
    print('''---------------------------------------------------------------------------------------------------------
    ******************************************** Set up the connections ********************************************
    ---------------------------------------------------------------------------------------------------------
    ''')
    # *****************************************
    #              Spike Threshold
    # *****************************************
    for ind, cellname in enumerate(population_names) :
        for n in range(network.populations[cellname].POP_SIZE) :
            network.populations[cellname].cells[n].create_spike_detector(**SpikeDetectionParam[ind])
    
    # *****************************************
    #              PP Connections
    # *****************************************
    PPconnectivity, synPPlist, PPspktrain = get_stim_Network(network, PPStimParameters, PPsynParameters,
                                                             PPStimNetworkConnParameters, deadcell=DeadCell) 
    
    # *****************************************
    #          MC Spontaneous Activity
    # *****************************************
    if SpontMCFactor :
        MCSpontActivity = MCSpontaneousActivity(network, DeadMCList=DeadMCList)
        
    # *****************************************
    #              Light Stimulation to inspect chaotic dynamics 
    # *****************************************
    if lightstim : 
        LightStimconnectivity, synLightStimlist, LightStimspktrain = get_stim_Network(network, LightStimParameters, LightStimsynParameters,
                                                                 LightStimNetworkConnParameters, deadcell=[DeadCell[0]])
    
    # *****************************************
    #              Cells Connection
    # *****************************************
    Networkconnectivity = [dict(celltype=name,connectivity=[]) for name in population_names]
    conncounttot = [[],[],[],[]]
    syncounttot = [[],[],[],[]]
    PathologicalNetwork = [dict(name='AMPAR',connectivity=[]),dict(name='KAR',connectivity=[])]
    # Delete the precedent network
    for i, pre in enumerate(population_names):
            for j, post in enumerate(population_names):
                print("------------ " + pre + "--->" + post + " ------------")
                # boolean connectivity matrix between pre- and post-synaptic
                # neurons in each population (postsynaptic on this RANK)
                npresyncell = network.populations[pre].POP_SIZE
                npostsyncell = network.populations[post].POP_SIZE
                connectivity = get_connectivity(
                    network,
                    pretype=pre, posttype=post,
                    npresyncell=npresyncell, npostsyncell=npostsyncell,
                    networkconnparam=networkConnParameters[i][j],
                    networktype=networktype,randomtype=randomtype,ringoption=ringoption
                )
                #Apply cell loss
                connectivity[DeadCell[i],:] = 0
                connectivity[:,DeadCell[j]] = 0
                Networkconnectivity[i]['connectivity'].append(connectivity)
                if pre == 'GC' and post == 'GC':
                    [indi,indj] = np.where(connectivity==True)
                    indindx = np.arange(0,len(indi))
                    KAindx = np.random.choice(indindx,int(KARfactor*len(indi)),replace=False)
                    AMPAindx = np.array([idx for idx in indindx if idx not in np.unique(KAindx)])
                    synindx = [AMPAindx,KAindx]
                    for typesyn in range(2) :
                        connectivitysyn = np.zeros(np.shape(connectivity))
                        for x in synindx[typesyn] :
                            connectivitysyn[indi[x]][indj[x]] = 1
                        connectivitysyn = connectivitysyn.astype(bool)
                        PathologicalNetwork[typesyn]['connectivity'] = connectivitysyn
                        (conncount, syncount) = network.connect(
                            pre=pre, post=post,
                            connectivity=connectivitysyn,
                            syntype=synapseModel,
                            synparams=GCGCsynParameters[typesyn],
                            synname=PathologicalNetwork[typesyn]['name'],
                            weightfun=weightFunction,
                            weightargs=GCGCweightArguments[typesyn],
                            minweight=minweight,
                            delayfun=delayFunction,
                            delayargs=delayArguments[i][j],
                            mindelay=mindelay,
                            multapsefun = None,
                            syn_pos_args=synapsePositionArguments[i][j],
                            save_connections=True,
                        )
                else :
                     (conncount, syncount) = network.connect(
                         pre=pre, post=post,
                         connectivity=connectivity,
                         syntype=synapseModel,
                         synparams=synapseParameters[i][j],
                         weightfun=weightFunction,
                         weightargs=weightArguments[i][j],
                         minweight=minweight,
                         delayfun=delayFunction,
                         delayargs=delayArguments[i][j],
                         mindelay=mindelay,
                         multapsefun = None,
                         syn_pos_args=synapsePositionArguments[i][j],
                         save_connections=True,
                     )
                conncounttot[i].append(conncount)
                syncounttot[i].append(syncount)
    
    
    print('''
    ---------------------------------------------------------------------------------------------------------
    *************************************** Measure options ***************************************
    ---------------------------------------------------------------------------------------------------------
    ''')
    # *****************************************
    #              LFP options
    #        Measure with 2 electrodes
    #               1 in the Hilus
    #                1 in the GCL
    # *****************************************
    numt = 100
    tlist = np.linspace(t0,t1,num = numt)
    tmid = minimize(tmaxresearch,t0+3, args=(t0,lengthGCL/2),tol=1e-6).x[0]
    #GCLline = np.array([r1*np.cos(tlist),r2*np.sin(tlist)]).T
    #centerpts = np.mean(GCLline,axis=0)
    centerGCL = np.array([r1*np.cos(tmid),r2*np.sin(tmid)])
    
    
    nbpts = 11 # Number of regular electrodes
    electrodeParametersLFP = dict(
        x=np.array([centerGCL[0]]*nbpts), 
        y=np.array([centerGCL[1]]*nbpts),
        z=np.linspace(0,6000,nbpts),
        n=50,  # nb of discrete point used to compute the potential
        sigma=.276,  # conductivity S/m
        method="linesource"
    )
    
    electrodeLFP = RecExtElectrode(cell=None, **electrodeParametersLFP)
    
    # *****************************************
    #  method Network.simulate() parameters
    # *****************************************
    networkSimulationArguments = dict(
        rec_pop_contributions=True,
        rec_imem=True,
        to_memory=False,
        to_file=True,
        file_name= 'LFP'+idnameIP
    )
    
    print('''
    ---------------------------------------------------------------------------------------------------------
    *************************************** Record the parameters used in this simulation ***************************************
    ---------------------------------------------------------------------------------------------------------
    ''')
    with open(os.path.join(OUTPUTPATH, 'Parameters'+idnameIP+'.txt'), "w") as f:
        f.write(' ***************************************** Experiment Options ***************************************** \n')
        for key,val in networkParameters.items() :
            f.write('{} = {}\n'.format(key,val))
        f.write('\n ***************************************** Network Options ***************************************** \n')
        for key,val in NetworkParameters2.items() :
            f.write('{} = {}\n'.format(key,val))
        f.write('\n ***************************************** PP Stimulation options ***************************************** \n')
        for key,val in PPStimParameters.items() :
            f.write('{} = {}\n'.format(key,val))
        f.write('\n ***************************************** PP Stimulation Parameters ***************************************** \n')
        f.write('\n --------- Synaptic Parameters --------- \n')
        for synparam in PPsynParameters :
            f.write('PostSynaptic')
            for key,val in synparam.items() :
                f.write('{} = {}\t'.format(key,val))
            f.write('\n')
        f.write('\n --------- Connectivity Parameters --------- \n')
        for namecell,param in zip(population_names,PPStimNetworkConnParameters) :
            f.write('######## {} ########\n'.format(namecell))
            for key,val in param.items() :
                f.write('\t{} = {}'.format(key,val))
            f.write('\n')
        f.write('\n ***************************************** Hilar Cell Loss ***************************************** \n')
        f.write('{}% of hilar cell loss \n'.format(int(100*hilarloss)))
        f.write('Cell loses : \n')
        for i, namecell in enumerate(population_names) :
            f.write('{} : \t'.format(namecell))
            for j in DeadCell[i] :
                f.write('{} \t'.format(j))
                f.write('\n')
        f.write('\n ***************************************** MC Spontaneous Activity ***************************************** \n')
        if SpontMCFactor :
            for nMC,MCSpont in enumerate(MCSpontActivity) :
                f.write('MC #{}'.format(nMC))
                for key,val in MCSpont.items() :
                    f.write('\t{} = {}'.format(key,val))
        else :
            f.write('No Spontaneous Activity of MCs \n')
        
        f.write('\n ***************************************** Cells models ***************************************** \n')
        for namecell,param in zip(population_names,CelltypeParameters) :
            f.write('######## {} ########\n'.format(namecell))
            for key,val in param.items() :
                f.write('\t{} = {}'.format(key,val))
            f.write('\n')
        f.write('\n ***************************************** Network Parameters ***************************************** \n')
        for precelltype,synparam,wparam,delayparam,synposparam,connparam in zip(population_names,
                                        synapseParameters,weightArguments,delayArguments,
                                        synapsePositionArguments, networkConnParameters) :
            f.write('######## Presyn : {} ########\n'.format(precelltype))
            for postcelltype,w,syn,delay,synpos,conn in zip(population_names,synparam,wparam,
                                        delayparam,
                                        synposparam,connparam) :
                f.write('\n ------------{}->{} ------------ \n'.format(precelltype,postcelltype))
                f.write('\t Synaptic Parameters')
                f.write('\n \t')
                for key,val in syn.items() :
                    f.write('\t{} = {}'.format(key,val))
                f.write('\n')
                f.write('\t Weight Parameters')
                f.write('\n \t')
                for key,val in w.items() :
                    f.write('\t{} = {}'.format(key,val))
                f.write('\n')
                f.write('\t Delay Parameters')
                f.write('\n \t')
                for key,val in delay.items() :
                    f.write('\t{} = {}'.format(key,val))
                f.write('\n')
                # f.write('\t Multapse Parameters')
                # f.write('\n \t')
                # for key,val in multaps.items() :
                #     f.write('\t{} = {}'.format(key,val))
                # f.write('\n')
                f.write('\t Synapse Position Parameters')
                f.write('\n \t')
                for key,val in synpos.items() :
                    f.write('\t{} = {}'.format(key,val))
                f.write('\n')
                f.write('\t Connectivity Parameters')
                f.write('\n \t')
                for key,val in conn.items() :
                    f.write('\t{} = {}'.format(key,val))    
                f.write('\n')
    


    print('''
    ---------------------------------------------------------------------------------------------------------
    *************************************** Set the stimulation trains ***************************************
    ---------------------------------------------------------------------------------------------------------
    ''')

    ActivePP = np.random.choice(np.arange(npp),size=(int(GCstim),PPnspk),replace=False)
    # Cancel all preexisting spike trains
    for synpp in synPPlist :
        for nsyn in range(len(synpp)) :
            for npost in range(len(synpp[nsyn])) :
                synpp[nsyn][npost].set_spike_times(np.zeros(0)) # np.zeros(0) = array([])
                
    PPspktrain1wind = [[] for ii in range(len(PPspktrain))]            
    for nwind in range(PPnspk) :
        for nsyn in ActivePP[:,nwind] :
            PPspktrain1wind[nsyn].append(PPspktrain[nsyn][nwind])
            
    PPspktrain1wind = [np.array(_) for _ in PPspktrain1wind]                         
    for synpp in synPPlist :
        for nwind in range(PPnspk) :
            for nsyn in ActivePP[:,nwind] :
                for npost in range(len(synpp[nsyn])) :
                    synpp[nsyn][npost].set_spike_times(PPspktrain1wind[nsyn])
                        
    if lightstim :                 
        ActiveLightStim = np.random.choice(np.arange(npp),size=(1),replace=False) # Activation of 1 random GC !!! A VERIFIER !!!
        for synlight in synLightStimlist :
           for nsyn in range(len(synlight)) :
               for npost in range(len(synlight[nsyn])) :
                   synlight[nsyn][npost].set_spike_times(np.zeros(0))      
                   if nsyn in ActiveLightStim :
                       synlight[nsyn][npost].set_spike_times(LightStimspktrain[nsyn])  
            

    print('''
    ---------------------------------------------------------------------------------------------------------
    *************************************** Start simulation and record the spikes ***************************************
    ---------------------------------------------------------------------------------------------------------
    ''')
    print(datetime.datetime.today().strftime("%Y-%m-%d %H:%M"))
    print(idname)
    
    print('\n --------- IP # {} --------- \n'.format(IP))
    
    if simvar :
        networkSimulationArguments['file_name']='LFP'+idnameIP
        
        SPIKES = network.simulate(
                probes=[electrodeLFP],
                **networkSimulationArguments
            )
        VcellArrays = [np.empty((network.populations[celltypepost].POP_SIZE, len(network.populations[celltypepost].cells[0].somav)),
                                dtype=[('vsoma', 'f8')]) for celltypepost in population_names]
        for ii,celltypepost in enumerate(population_names) :
            for jj in range(network.populations[celltypepost].POP_SIZE) :
                VcellArrays[ii][jj]['vsoma'] = network.populations[celltypepost].cells[jj].somav
        with h5py.File(os.path.join(OUTPUTPATH2,'VcellArrays'+idnameIP+'.h5') ,'a') as f:
            for ii,celltypepost in enumerate(population_names) :
                if celltypepost in f.keys():
                    del f[celltypepost]
                f[celltypepost] = VcellArrays[ii]
                
        with open(os.path.join(OUTPUTPATH2, 'SPIKES'+idnameIP), "wb") as f:
            pickle.dump(SPIKES,f)
            
        
        with open(os.path.join(OUTPUTPATH2, 'OUTPUT'+idnameIP+'.txt'), 'w') as f:
            for jj in range(len(population_names)) : 
                for spt, gid in zip(SPIKES['times'][jj], SPIKES['gids'][jj]):
                    for ii in range(len(spt)) :
                        f.write('{}\t{}'.format(spt[ii],gid))
                        f.write('\n')
            
        PPspktraintemp = [PPspktrain1wind[ii] for ii in np.unique(ActivePP)]
        with open(os.path.join(OUTPUTPATH2, 'INPUT'+idnameIP+'.txt'), 'w') as f:
            for indpp,ppspktrain in zip(np.unique(ActivePP),PPspktraintemp) :
                for ppspk in ppspktrain :
                    f.write('{}\t{}'.format(ppspk,indpp))
                    f.write('\n')
        if lightstim : 
            with open(os.path.join(OUTPUTPATH2, 'LightINPUT'+idnameIP+'.txt'), 'w') as f:
                for indstim,stimspktrain in zip(np.unique(ActiveLightStim),LightStimspktrain) :
                    for stimspk in stimspktrain :
                        f.write('{}\t{}'.format(ppspk,indstim))
                        f.write('\n')
        f.close()    
    
        
        '''
        ----------------------------------------------------------------------------------------------------------------------------
        ****************************************************** Display results ******************************************************
        ----------------------------------------------------------------------------------------------------------------------------
        '''
        
        colorsme = dict(GC=[.25,.45,.9],
                  MC=[.95,.6,.15],
                  BC=[.75,.25,.25],
                  HC=[.15,.75,.25],
                  PP=[153/256,0,153/256])
        
        
        def remove_axis_junk(ax, lines=['right', 'top']):
            """remove chosen lines from plotting axis"""
            for loc, spine in ax.spines.items():
                if loc in lines:
                    spine.set_color('none')
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
        
        
        def draw_lineplot(
                ax, data, dt=0.1,
                T=(0, 200),
                scaling_factor=1.,
                vlimround=None,
                label='local',
                scalebar=True,
                scalebarunit=False,
                unit='mV',
                ylabels=True,
                xlabels=True,
                remaxisjunk=True,
                color='r',
                ztransform=True,
                filter_data=False,
                filterargs=dict(N=2, Wn=0.02, btype='lowpass')):
            """helper function to draw line plots"""
            tvec = np.arange(data.shape[1]) * dt
            tinds = (tvec >= T[0]) & (tvec <= T[1])
        
            # apply temporal filter
            if filter_data:
                b, a = ss.butter(**filterargs)
                data = ss.filtfilt(b, a, data, axis=-1)
        
            # subtract mean in each channel
            if ztransform:
                dataT = data.T - data.mean(axis=1)
                data = dataT.T
        
            zvec = -np.arange(data.shape[0])
            vlim = abs(data[:, tinds]).max()
            if vlimround is None:
                vlimround = 10.**np.round(np.log(vlim)/np.log(10)) / scaling_factor
            else:
                pass
            yticklabels = []
            yticks = []
        
            for i, z in enumerate(zvec):
                if i == 0:
                    ax.plot(tvec[tinds], data[i][tinds] / vlimround + z, lw=0.75,
                            rasterized=False, label=label, clip_on=False,
                            color=color)
                else:
                    ax.plot(tvec[tinds], data[i][tinds] / vlimround + z, lw=0.75,
                            rasterized=False, clip_on=False,
                            color=color)
                yticklabels.append('ch. %i' % (i + 1))
                yticks.append(z)
        
            if scalebar:
                ax.plot([tvec[-1], tvec[-1]],
                        [-1, -2], lw=2, color='k', clip_on=False)
                if scalebarunit :
                    ax.text(tvec[-1] + np.diff(T) * 0.02, -1.5,
                            f'{vlimround}'+ f'{unit}',
                            color='k', rotation='vertical',
                            va='center')
                else :
                    powerten = int(np.floor(np.log10(vlimround)))
                    ax.text(tvec[-1] + np.diff(T) * 0.02, -1.5,
                        '$'+f'{int(vlimround/10**powerten)}'+'.10^{' + f'{powerten}' + '}$ ' + f'{unit}',
                        color='k', rotation='vertical',
                        va='center')
        
            ax.axis(ax.axis('tight'))
            ax.yaxis.set_ticks(yticks)
            if ylabels:
                ax.yaxis.set_ticklabels(yticklabels)
                ax.set_ylabel('channel', labelpad=0.1)
            else:
                ax.yaxis.set_ticklabels([])
            if remaxisjunk :
                remove_axis_junk(ax, lines=['right', 'top'])
            if xlabels :
                ax.set_xlabel(r't (ms)', labelpad=0.1)
        
            return vlimround
        
        
        # ********************************************
        #                 PP activity
        # ********************************************
        if nspr == 0 :
            titlebase = 'Healthy Network'
        else :
            if KARfactor > 0 :
                titlebase = '{}%A/{}%K%' .format(int((1-KARfactor)*100),int(KARfactor*100))
            else :
                titlebase = '100%AMPAR'
            titlebase+= ' - {}%spr - {}%PPint - {}GCs stim'.format(nspr,PPintervals,int(GCstim*npp))
            titlebase+='{}'.format(['Base','LightStim'][lightstim>0])
            
        f,ax = plt.subplots(dpi=200)
        yticks = []
        PPspktraintemp = [PPspktrain1wind[ii] for ii in np.unique(ActivePP)]
        for nPP,PPspk in enumerate(PPspktraintemp):
            for t in PPspk :
                ax.plot(t,nPP,'o',color=colorsme['PP'],ms=3)
            yticks.append(nPP)
            
        remove_axis_junk(ax, lines=['right', 'top'])
        ax.set_title('PP activity'+' - '+titlebase)
        ax.yaxis.set_ticks(yticks)
        ax.set_xlabel('t (ms)')
        ax.set_ylabel('#PPTrain')
        ax.set_xlim([0,tstop+tstop//100])
        # if tstop<1000 :
        #     ax.xaxis.set_ticks(np.arange(0,tstop+100,100))
            # ax.set_xticklabels(['' for i in range(tstop//100+1)])
        if saveopt :
            plt.savefig(join(OUTPUTPATH,'PPactivity-'+idnameIP+'.png'), dpi=200, bbox_inches='tight')
            # plt.savefig(join(OUTPUTPATH,'PPactivity-'+idnameIP+imgfrm), dpi=200, bbox_inches='tight')
            plt.close(f)
        
        
        # ********************************************
        #            Activity Rasterplot
        # ********************************************
        
        fig, ax = plt.subplots(2, 1,gridspec_kw={'height_ratios': [9,1]},dpi=200)
        for name, spts, gids in zip(
                population_names, SPIKES['times'], SPIKES['gids']):
            t = []
            g = []
            for spt, gid in zip(spts, gids):
                t = np.r_[t, spt]
                g = np.r_[g, np.zeros(spt.size) + gid]
            ax[0].plot(t, g, 'o',color=colorsme[name], ms=1, label=name)
            # ax.plot(t[t >= 200], g[t >= 200], '.',color=colorsme[name], ms=3, label=name)
        ax[0].legend(loc=1)
        remove_axis_junk(ax[0], lines=['right', 'top'])
        ax[0].set_ylabel('#Cell')
        ax[0].set_title('Rasterplot'+' - '+titlebase)
        ax[0].set_xlim([0,tstop+tstop//10])
        # if tstop<1000 :
        #     ax[0].xaxis.set_ticks(np.arange(0,tstop+100,100))
        #     ax[0].set_xticklabels(['' for i in range(tstop//100+1)])
        
        for nPP,PPspk in enumerate(PPspktraintemp):
            for t in PPspk :
                ax[-1].plot(t,nPP,'o',color=colorsme['PP'],ms=2)
        # if tstop<1000 :
        #     ax[-1].xaxis.set_ticks(np.arange(0,tstop+100,100))
        ax[-1].set_xlim([0,tstop+tstop//10])
        ax[-1].yaxis.set_ticks([])
        ax[-1].set_ylim(np.array([-1,len(np.unique(ActivePP))]))
        ax[-1].set_xlabel('t (ms)')
        ax[-1].set_ylabel('PP')
        if saveopt :
            plt.savefig(join(OUTPUTPATH,'Rasterplot-'+idnameIP+'.png'), dpi=200, bbox_inches='tight')
            # plt.savefig(join(OUTPUTPATH,'Rasterplot-'+idnameIP+imgfrm), dpi=200, bbox_inches='tight')
            plt.close(fig)
        
        # ********************************************
        #                  LFP
        # ********************************************
        
        filterLFPvar = 8 #16 in LFPy examples
        lfpname = os.path.join(OUTPUTPATH2, 'LFP'+idnameIP)
        with h5py.File(lfpname, 'r') as f:
            print(f['RecExtElectrode0'])
            lfpalltype = f['RecExtElectrode0'][population_names[0]]
            for name in population_names[1:] :
                lfpalltype += f['RecExtElectrode0'][name]                
            lfp = ss.decimate(lfpalltype, q=filterLFPvar,
                              zero_phase=True)
            
        
        fig, ax = plt.subplots(1, 1,dpi=200)      
        draw_lineplot(ax,
                      lfp,
                      dt=network.dt * filterLFPvar,
                      T=(0, tstop+tstop//10),
                      scaling_factor=.5,
                      vlimround=2e-3,
                      label='LFP',
                      scalebar=True,
                      unit='mV',
                      ylabels=True,
                      xlabels=True,
                      color=[0,0,0],#f'C{i}',
                      ztransform=True
                      )
        ax.set_title('LFP in GCL'+' - '+titlebase)
        ax.set_xlim([0,tstop+tstop//10])
        # if tstop<1000 :
        #     ax.xaxis.set_ticks(np.arange(0,tstop+100,100))
        #     ax.set_xticklabels(['' for i in range(tstop//100+1)])
        
        if saveopt :
            plt.savefig(join(OUTPUTPATH,'LFPGCL-'+idnameIP+'.png'), dpi=200, bbox_inches='tight')
            # plt.savefig(join(OUTPUTPATH,'LFPGCL-'+idnameIP+imgfrm), dpi=200, bbox_inches='tight')
            plt.close(fig)
        # ax.set_xlim([0,tstop*5])
        
        ''' ---Single Electrode LFP --- '''                
        fig, ax = plt.subplots(1, 1,dpi=200)      
        draw_lineplot(ax,
                      np.array([lfp[len(lfp)//2,:],np.nan*np.empty(len(lfp[0]))]),
                      dt=network.dt * filterLFPvar,
                      T=(0, tstop+tstop//10),
                      scaling_factor=.5,
                      vlimround=2e-3,
                      label='LFP',
                      scalebar=True,
                      unit='mV',
                      ylabels=True,
                      xlabels=True,
                      color=[0,0,0],#f'C{i}',
                      ztransform=True
                      )
        ax.set_title('LFP in GCL'+' - '+titlebase)
        ax.set_xlim([0,tstop+tstop//10])
        # if tstop<1000 :
        #     ax.xaxis.set_ticks(np.arange(0,tstop+100,100))
        #     ax.set_xticklabels(['' for i in range(tstop//100+1)])
        
        if saveopt :
            plt.savefig(join(OUTPUTPATH,'LFPGCL-'+idnameIP+'.png'), dpi=200, bbox_inches='tight')
            plt.savefig(join(OUTPUTPATH,'LFPGCL-'+idnameIP+imgfrm), dpi=200, bbox_inches='tight')
            # plt.savefig(join(OUTPUTPATHValerie,'LFPGCL-'+idnameIP+'.png'), dpi=200, bbox_inches='tight')
            # plt.savefig(join(OUTPUTPATHValerie,'LFPGCL-'+idnameIP+imgfrm), dpi=200, bbox_inches='tight')
            plt.close(fig)
        # ax.set_xlim([0,tstop*5])
        
        
        # ********************************************
        #         Soma potentials of some cells
        # ********************************************
        ncellpick = [np.random.randint(0,ngcell,25), np.random.randint(0,nmcell,3), np.random.randint(0,nbcell,2), np.random.randint(0,nhcell,2)]
        scalebarvar = [True,False,False,False]
        t = np.linspace(0,len(network.populations['GC'].cells[0].somav)*dt,len(network.populations['GC'].cells[0].somav))
        
        fig, ax = plt.subplots(nrows=5, ncols=1,gridspec_kw={'height_ratios': [7.35,.7,.7,.7,.7]},dpi=200,sharex=False)
        
        for i, cellname in enumerate(population_names) :
            somavs_pop = np.array([network.populations[cellname].cells[n].somav for n in ncellpick[i]])  # avoid undeclared variable            
            draw_lineplot(ax[i],
                            somavs_pop,
                            dt=dt,
                            T=(0, tstop+tstop//10),
                            scaling_factor=.25,
                            vlimround=70,
                            label=cellname,
                            scalebar=scalebarvar[i],
                            scalebarunit = 70,
                            unit='mV',
                            ylabels=True,
                            xlabels=True,
                            remaxisjunk=True,
                            color=colorsme[cellname],
                            ztransform=True
                            )
            ax[i].set_ylabel(cellname)
            ax[i].yaxis.set_ticks([])
            ax[i].set_xlabel('')
            ax[i].set_xlim([0,tstop+tstop//10])
            ax[i].set_xticklabels('')
            # if tstop<1000 :
            #     ax[i].xaxis.set_ticks(np.arange(0,tstop+100,100))
            #     ax[i].set_xticklabels(['' for i in range(tstop//100+1)])
        ax[0].set_title('Soma potential'+' - '+titlebase)
        # if tstop<1000 :
        #     ax[0].xaxis.set_ticks(np.arange(0,tstop+100,100))
        #     ax[0].set_xticklabels(['' for i in range(tstop//100+1)])
        
        for nPP,PPspk in enumerate(PPspktraintemp):
            for t in PPspk :
                ax[-1].plot(t,nPP,'o',color=colorsme['PP'],ms=1.5)
        # if tstop<1000 :
        #     ax[-1].xaxis.set_ticks(np.arange(0,tstop+100,100))
        ax[-1].set_xlim([0,tstop+tstop//10])
        ax[-1].yaxis.set_ticks([])
        ax[-1].set_ylim([-1,len(np.unique(ActivePP))])
        ax[-1].set_xlabel('t (ms)')
        ax[-1].set_ylabel('PP')
        # remove_axis_junk(ax[-1], lines=['right', 'top'])
        if saveopt :
            plt.savefig(join(OUTPUTPATH,'CellsSomaPotential-'+idnameIP+'.png'), dpi=200, bbox_inches='tight')
            # plt.savefig(join(OUTPUTPATH,'CellsSomaPotential-'+idnameIP+imgfrm), dpi=200, bbox_inches='tight')
            plt.close(fig)
        
        
        ax = plt.figure(dpi=dpi,figsize=(10,12)).add_subplot(projection='3d')
        idx = 0
        for k, cellname in enumerate(population_names) :
            xcell = []
            ycell = []
            zcell = []
            sz = []
            for indcell in range(population_sizes[k]) :
                if cellname == 'GC' and ((indcell>125) or (indcell<85))  :
                    continue
                xcell.append(np.mean(network.populations[cellname].cells[indcell].x[idx]))
                ycell.append(np.mean(network.populations[cellname].cells[indcell].y[idx]))
                zcell.append(np.mean(network.populations[cellname].cells[indcell].z[idx]))
                sz = 1*network.populations[cellname].cells[0].d[idx]
            ax.scatter(xcell,ycell,zcell,s=sz,marker='o',
                       c=colorsme[cellname], label=cellname)

        # Draw electrode
        ax.plot(electrodeParametersLFP['x'],electrodeParametersLFP['y'],electrodeParametersLFP['z'],'o',
                color=[255/255, 153/255, 204/255],ms=3.5,label='LFP Electrode')
        ax.set_xlabel(r'X($\mu$m)')
        ax.set_ylabel(r'Y($\mu$m)')
        ax.set_zlabel(r'Z($\mu$m)')
        ax.legend(loc=3)
        ax.set_title('3D Geometry of some cells'+' - '+titlebase)
        
        #zlimtemp = np.ceil(np.array(ax.get_zlim())*10)//10
        ax.set(xlim=(-1000, -100))
        ax.set(ylim=(-150,300))
        ax.set(zlim=(1000,1500))
        # ax.set(zlim=(zlimtemp[0], zlimtemp[1]))
        ax.grid(False)
        
        elev = 10
        azim = -60
        roll = 0
        ax.view_init(elev=elev, azim=azim, roll=roll)
        if saveopt :
            plt.savefig(join(OUTPUTPATH,'NetworkGeometrySomaOnlyWithAxis-'+idname+'.png'), dpi=200, bbox_inches='tight')
        # Scale 
        ax.plot([-250,-150],[0,0],[0,0],'k',lw=2) # x-scale
        ax.plot([-150,-150],[-100,0],[0,0],'k',lw=2) # y-scale
        ax.plot([-150,-150],[0,0],[0,1000],'k',lw=2) # z-scale
        ax.text(-300,0,0,'x')
        ax.text(-150,-150,0,'y')
        ax.text(-150,0,1100,'z')
        
        scaletxt = r'$ x : 100 \mu m$'+' \n' +r'$ y : 100 \mu m $'+' \n' +r'$ z : 1000 \mu m $'
        ax.text(0,-100,-500,scaletxt)
        
        ax.set_axis_off()
        ax.set_box_aspect(None, zoom=1.2)


        if saveopt :
            plt.savefig(join(OUTPUTPATH,'NetworkGeometrySomaOnlyWithoutAxis-'+idname+'.png'), dpi=200, bbox_inches='tight')
            
        # ********************************************
        #          3D Geometry of some cells
        # ********************************************
        
        ax = plt.figure(dpi=200).add_subplot(projection='3d')
        for k, cellname in enumerate(population_names) :
            for indcell in ncellpick[k] :
                for sec in network.populations[cellname].cells[indcell].allseclist:
                    idx = network.populations[cellname].cells[indcell].get_idx(sec.name())
                    for i in idx:
                        ax.plot(network.populations[cellname].cells[indcell].x[i], 
                                network.populations[cellname].cells[indcell].y[i], 
                                network.populations[cellname].cells[indcell].z[i],
                                color=colorsme[cellname],
                                lw=network.populations[cellname].cells[indcell].d[i],)
        # Draw electrode
        ax.plot(electrodeParametersLFP['x'],electrodeParametersLFP['y'],electrodeParametersLFP['z'],'o',
                color=[255/255, 153/255, 204/255],ms=3.5,label='LFP Electrode')
        ax.set_xlabel(r'X($\mu$m)')
        ax.set_ylabel(r'Y($\mu$m)')
        ax.set_zlabel(r'Z($\mu$m)')
        ax.legend()
        ax.set_title('3D Geometry of some cells'+' - '+titlebase)
        if saveopt :
            plt.savefig(join(OUTPUTPATH,'SomeNetworkCellsGeometry-'+idname+'.png'), dpi=200, bbox_inches='tight')
        
        '''
        ---------------------------------------------------------------------------------------------------------
        ******************************************** End of experiment ********************************************
        ---------------------------------------------------------------------------------------------------------
        '''
        
    network.pc.gid_clear()
    return '--------- '+datetime.datetime.today().strftime("%Y-%m-%d %H:%M")+'\n'+str((IP,nspr,PPintervals,KARfactor,PPnspk,resultsfolderpathname,GCstim))+' DONE ---------'
    
    
    
# !!! Before using the code you have to select the right folders if you don't want in the same folder !
cellpath = os.getcwd() #  Path of the folder where your neuron models are
scriptpath = os.getcwd() #  Path of the folder where your script is
resultsfolderpathname = cellpath # Path of the folder where you want to store the results

ParametersSimulation = dict()  
with open(os.path.join(cellpath, 'ParametersSimulationSingle'+'.txt'), "r") as f:  
    for line in f:
        temp = line.split('\t')
        typval = temp[2]
        key = temp[0]
        if "<class 'int'>" in typval:   
            ParametersSimulation[key]=int(temp[1])
        elif "<class 'bool'>" in typval:
            ParametersSimulation[key]=bool(temp[1])
        elif "<class 'float'>" in typval:
            ParametersSimulation[key]=float(temp[1])
        elif "<class 'list'>" in typval:
            ParametersSimulation[key]=[int(x) for x in temp[1][1::3]]
    #trial = ParametersSimulation['trial']


try :
    neuron.load_mechanisms(cellpath)
except :
    pass

os.chdir(cellpath)
resultsfolderpathname =  cellpath 
    

    
        
#!!! Lauch the simulation
Network500PPs_1IPfor1window_1PP1GC_PPintervalExploration(ParametersSimulation['IP'],ParametersSimulation['nspr'],
                                                         ParametersSimulation['PPintervals'], 0,ParametersSimulation['PPnspk'],
                                                         resultsfolderpathname, cellpath, ParametersSimulation['GCstim'], 
                                                         ParametersSimulation['nGClayer'], ParametersSimulation['simvar'],
                                                         ParametersSimulation['lightstim'],ParametersSimulation['randseed'])                                                         