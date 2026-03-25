TITLE Ihfast channel TITLE Ihfast.mod:
:Ihfast channel (HCN hyperpolarization-activated [cyclic nucleotide-gated] cation channel)
 
COMMENT
Fast part of the Ih 
Ihfast.mod adapted from  :
hyperde3.mod (hyf)
Modified from Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
https://modeldb.science/155568?tab=2&file=TejadaEtAl2014/hyperde3.mod
and HCN.mod (hyf)
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
SUFFIX ihfast 
NONSPECIFIC_CURRENT ihf
RANGE ghf
RANGE ghfbar, ehf, ihf
RANGE hfinf, hftau
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
	
	ehf = -40 (mV)
	
	ghfbar = 0 (mho/cm2)	
}


ASSIGNED {
	ghf (mho/cm2)
  
	ihf (mA/cm2)

	hfinf 
 	hftau (ms)
	
	tadj
} 
 
STATE {
	hf 
}


INITIAL {
	rates(v)
	
	hf = hfinf	
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ghf = ghfbar * hf*hf	
	ihf = ghf * (v-	ehf)
		}
 

DERIVATIVE states {	
	rates(v)	:      at the current v and dt.
	
	hf' = (hfinf-hf)/hftau
}
 

PROCEDURE rates(vm (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	tadj = q10^((celsius - temp)/10)
       
	:"hf" FAST CONTROL Hype activation system
	hfinf =  1 / (1 + exp( (vm+91)/10 ))
	hftau = (14.9 + 14.1 / (1+exp(-(vm+95.2)/0.5)))/tadj 
}
 
 

