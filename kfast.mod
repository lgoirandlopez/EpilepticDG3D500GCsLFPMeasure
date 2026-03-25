TITLE kfast.mod
: Sodium channel, Hodgkin-Huxley style kinetics.  

COMMENT
Original Mod Files:
Original names 'ichan2.mod' 
Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
https://modeldb.science/155568?tab=2&file=TejadaEtAl2014/ichan2.mod
Original modified from Santhakumar 2005 

Author: Lucas Goirand-Lopez
Inspired from Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kfast
	USEION k READ ek WRITE ik VALENCE 1
	RANGE nf, gkf, gkfbar
	RANGE nfinf, nftau
	RANGE vshift, vcoef
	GLOBAL q10, temp, tadj, vmin, vmax
}

PARAMETER {
	gkfbar = 0.016   	(S/cm2)	: 0.12*1e4 (pS/um2)
								

	temp = 6.3	(degC)		: original temp 
	q10  = 3				: temperature sensitivity

	v 		(mV)
	celsius	= 6.3	(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
	
	vshift = 0 (mV)
	vcoef = 1 (mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 		(mA/cm2)
	gkf		(S/cm2)
	ek		(mV)
	nfinf
	nftau 	(ms)
	tadj
}
 

STATE { nf }

INITIAL { 
	rates(v)
	nf = nfinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gkf = gkfbar*nf*nf*nf*nf
	ik = gkf * (v - ek)
} 


DERIVATIVE states {    
	rates(v)
	nf' = (nfinf - nf)/nftau
}



PROCEDURE rates(vm (mV)) {  

	:TABLE  nftau, nfinf DEPEND celsius FROM vmin TO vmax WITH 199
	
	LOCAL  alpha, beta, sum

	alpha = -0.07*vtrap((vm+65-47+vshift),-6)
	beta = 0.264/exp((vm+65-22+vshift)/40)
	tadj = q10^((celsius - temp)/10)
	nftau = 1/(alpha+beta)/tadj      
	nfinf = alpha/(alpha+beta)	
}


FUNCTION vtrap(x (mV), y (mV)) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}



