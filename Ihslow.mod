TITLE Ihslow.mod
:Ihslow channel (HCN hyperpolarization-activated [cyclic nucleotide-gated] cation channel)
 
COMMENT
Slow part of the Ih 
Ihfast.mod adapted from  :
hyperde3.mod (hys)
Modified from Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
https://modeldb.science/155568?tab=2&file=TejadaEtAl2014/hyperde3.mod
and HCN.mod (hys)
Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
https://modeldb.science/185355?tab=2&file=Yim_et_al_2015/HCN.mod

Original Mod File:
Original name 'hyperde3.mod'
Santhakumar V, Aradi I, Soltesz I (2005) J Neurophsiol 93:437-53 
https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=51781&file=/dentategyrusnet2005/hyperde3.mod

Mod File history:
Chen K, Aradi I, Thon N, Eghbal-Ahmadi M, Baram TZ, Soltesz I: Persistently
modified h-channels after complex febrile seizures convert the seizure-induced
enhancement of inhibition to hyperexcitability. Nature Medicine, 7(3) pp. 331-337, 2001.
(modeling by Ildiko Aradi, iaradi@uci.edu)
distal dendritic Ih channel kinetics for both HT and Control animals

Author: Lucas Goirand-Lopez 
Inspired from Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
And Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
ENDCOMMENT
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)} 
 
NEURON { 
SUFFIX ihslow 
NONSPECIFIC_CURRENT ihs
RANGE ghs, ehs, ihs
RANGE ghsbar
RANGE hsinf, hstau
GLOBAL q10, temp, tadj, vmin, vmax
}
 
UNITS {
        (mA) =(milliamp)
        (mV) =(millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
	
}
 
 
PARAMETER {
	v (mV) 
	temp = 6.3	(degC)		: original temp 
	q10  = 3				: temperature sensitivity
	celsius (degC)
	
	vmin = -120	(mV)
	vmax = 100	(mV)
	
	ehs = -40 (mV)
	
	ghsbar = 0 (mho/cm2)
}


ASSIGNED {
 	ghs (mho/cm2)
	
	ihs (mA/cm2)
	
	hsinf
 	hstau (ms)
	
	tadj
} 
 
STATE {
	hs
}


INITIAL {
	rates(v)
	
	hs = hsinf	
}

BREAKPOINT {	
	SOLVE states METHOD cnexp
	ghs = ghsbar * hs*hs
	ihs = ghs * (v-ehs)	
		}
 

DERIVATIVE states {	
	rates(v)	:      at the current v and dt.
	
	hs' = (hsinf-hs)/hstau
}
 

PROCEDURE rates(vm (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	tadj = q10^((celsius - temp)/10)

	:"hs" SLOW CONTROL Hype activation system
	hsinf =  1 / (1 + exp( (vm+91)/10 ))
	hstau = (80 + 172.7 / (1+exp(-(vm+59.3)/-0.83)))/tadj 
}
 
 

