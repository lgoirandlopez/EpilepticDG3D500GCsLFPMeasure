TITLE kir.mod
:Inward fast delayed rectifier potassium channel (Kir), Hodgkin-Huxley style kinetics

COMMENT
Original Mod Files:
Original names 'ichan2.mod' 
Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
https://modeldb.science/185355?tab=2&file=Yim_et_al_2015/Kir.mod

Mod File history:
- tau(V), linf(V) fitted to experimental values of human dentate gyrus granual cells
- ModelDB file adapted from 
  Wolf JA, Moyer JT, Lazarewicz MT, Contreras D, Benoit-Marand M, O'Donnell P, Finkel LH (2005) J Neurosci 25:9080-95
  https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=112834&file=/nacb_msp/kir.mod
- file modified to uses nomoclature of 
  Li X, Ascoli GA (2006) J of Comput Neurosci 21(2):191-209 
  Li X, Ascoli GA (2008) Neural Comput 20:1717-31

Author: Lucas Goirand-Lopez adapted from Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kir
	USEION k READ ek WRITE ik VALENCE 1
	RANGE  gkkirbar, vhalfl, kl, vhalft, at, bt, gkkir, l
	GLOBAL q10, temp, vmin, vmax
}

PARAMETER {
	v 		(mV)
	gkkirbar  = 1.44e-05	(S/cm2) 	: to be fitted     	

	: Boltzman steady state curve	
	vhalfl = -98.923594  (mV)    		: fitted to patch data, Stegen et al. 2012
	kl = 10.888538      (mV)    		: Stegen et al. 2012

	: tau_infty 
	vhalft=67.0828	 (mV)    		: fitted #100 \muM sens curr 350a,  Stegen et al. 2012
	at=0.00610779(/ms)   		: Stegen et al. 2012
	bt=0.0817741	 (/ms)	 		: Note: typo in Stegen et al. 2012

	: Temperature dependence
	temp = 6.3	(degC)		: original temp 
	q10  = 3				: temperature sensitivity

	celsius	= 6.3	(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
	
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 		(mA/cm2)
	gkkir		(S/cm2)
	ek		(mV)
	linf
	taul 	(ms)
}
 

STATE { l }

INITIAL { 
	rate(v)
	l=linf
}

BREAKPOINT {
	SOLVE states METHOD cnexp	: solve differential equations in states with method 'cnexp'
	gkkir = gkkirbar*l			: use state l to calulate gkkir
	ik = gkkir * ( v - ek )		: calculate ik 
} 


DERIVATIVE states {    
	rate(v)
	l' =  (linf - l)/taul		: differential equation 
}



PROCEDURE rate(vm (mV)) {  

	:TABLE  nftau, nfinf DEPEND celsius FROM vmin TO vmax WITH 199
	
	LOCAL tadj
	tadj=q10^((celsius-33)/10) 	
	linf = 1/(1 + exp((vm-vhalfl)/kl))			: l_steadystate
 	taul = 1/(tadj *(at*exp(-vm/vhalft) + bt*exp(vm/vhalft) ))	
}



