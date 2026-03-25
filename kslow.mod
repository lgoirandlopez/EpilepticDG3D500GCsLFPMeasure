TITLE kslow.mod
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
	SUFFIX kslow
	USEION k READ ek WRITE ik VALENCE 1
	RANGE ns, gks, gksbar
	RANGE nsinf, nstau
	RANGE vshift, vcoef
	GLOBAL q10, temp, tadj, vmin, vmax
}

PARAMETER {
	gksbar = 0.016   	(S/cm2)	: 0.12*1e4 (pS/um2)
								
	temp = 6.3	(degC)		: original temp 
	q10  = 3				: temperature sensitivity

	v 		(mV)
	celsius = 6.3		(degC)
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
	gks		(S/cm2)
	ek		(mV)
	nsinf
	nstau 	(ms)
	tadj
}
 

STATE { ns }

INITIAL { 
	rates(v)
	ns = nsinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gks = gksbar*ns*ns*ns*ns
	ik = gks * (v - ek)
} 


DERIVATIVE states {    
	rates(v)
	ns' = (nsinf - ns)/nstau
}



PROCEDURE rates(vm (mV)) {  

	:TABLE  nstau, nsinf DEPEND celsius FROM vmin TO vmax WITH 199
	
	LOCAL  alpha, beta, sum

	alpha = -0.028*vtrap((vm+65-35+vshift),-6)
	beta = 0.1056/exp((vm+65-10+vshift)/40)
	tadj = q10^((celsius - temp)/10)
	nstau = 1/(alpha+beta)/tadj      
	nsinf = alpha/(alpha+beta)	
}


FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}



