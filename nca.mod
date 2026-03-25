TITLE nca.mod 
: N-type calcium channel 

COMMENT
Original Mod Files:
Original names 'nca.mod' 
Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
https://modeldb.science/155568?tab=2&file=TejadaEtAl2014/nca.mod
Original modified from Migliore CA3  

Author: Lucas Goirand-Lopez
Adapted from Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
ENDCOMMENT
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX nca
	USEION nca READ enca WRITE inca VALENCE 2 
	RANGE gncabar, cinf, dinf , ctau, dtau, inca
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (milli/liter)

	FARADAY = (faraday) (coulomb)
	R = 8.313424  (joule/degC)
}

PARAMETER {
	v (mV)
	celsius	(degC)
	gncabar = .002 (mho/cm2)
	cai (mM)
	cao (mM)
	temp = 6.3	(degC)
	q10  = 3				: temperature sensitivity
}


STATE {
	c d 
}

ASSIGNED {
	inca (mA/cm2)
	:gnca (mho/cm2)
	cinf 		dinf
	ctau (ms)	dtau (ms)
	tadj
	enca (mV)
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	inca = gncabar*c*c*d*(v-enca)
}

UNITSOFF

DERIVATIVE states {	: exact when v held constant
	crate(v)
	drate(v)
	c' = (cinf - c)/ctau
	d' = (dinf - d)/dtau
}


INITIAL {
	crate(v)
	drate(v)
	c = cinf
	d = dinf
}


PROCEDURE crate(vm(mV)) {
	LOCAL  alpha, beta
	
	alpha = -0.19*vtrap(vm-19.88,-10)
	beta = 0.046*exp(-vm/20.73)
	tadj = q10^((celsius - temp)/10)
	
	cinf = alpha/(alpha+beta)
	ctau = 1/(alpha+beta)/tadj
}

PROCEDURE drate(vm(mV)) {
:
	LOCAL  alpha, beta
		
	alpha = 0.00016/exp(-vm/48.4)
	beta = 1/(exp((-vm+39)/10)+1)
	tadj = q10^((celsius - temp)/10)
	
	dinf = alpha/(alpha+beta)
	dtau = 1/(alpha+beta)/tadj
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON