TITLE T-calcium channel 
: T-type calcium channel

COMMENT
Original Mod Files:
Original names 'tca.mod' 
Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
https://modeldb.science/155568?tab=2&file=TejadaEtAl2014/tca.mod
Original modified from Migliore CA3  

Author: Lucas Goirand-Lopez
Adapted from Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX tca
	USEION tca READ etca WRITE itca VALENCE 2
	USEION ca READ cai, cao VALENCE 2
	RANGE gtcabar, ainf, binf , atau, btau
	RANGE itca, etca
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
	celsius = 6.3	(degC)
	gtcabar (mho/cm2)
	cai (mM)
	cao (mM)
}


STATE {
	a b 
}

ASSIGNED {
	itca (mA/cm2)
	ainf 		binf
	atau (ms)	btau (ms)
	etca (mV)
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	itca = gtcabar*a*a*b * ghk(v, cai, cao)
}

DERIVATIVE states {	: exact when v held constant
	arate(v)
	brate(v)
	a' = (ainf - a)/atau
	b' = (binf - b)/btau
}


INITIAL {
	arate(v)
	brate(v)
	a = ainf
	b = binf
	VERBATIM
	cai=_ion_cai;
	ENDVERBATIM
}


PROCEDURE arate(vm(mV)) {
	LOCAL  alpha, beta
	
	alpha = 0.2*(-1.0*vm+19.26)/(exp((-1.0*vm+19.26)/10.0)-1.0)
	beta = 0.009*exp(-vm/22.03)
	
	ainf = alpha/(alpha+beta)
	atau = 1/(alpha+beta)
}

PROCEDURE brate(vm(mV)) {
:
	LOCAL  alpha, beta
		
	alpha = 1.e-6*exp(-vm/16.26)
	beta = 1/(exp((-vm+29.79)/10.)+1.)
	
	binf = alpha/(alpha+beta)
	btau = 1/(alpha+beta)
}

FUNCTION ghk(vm(mV), ci(mM), co(mM)) (mV) {
	LOCAL nu,f

	f = (1000)*R*(celsius +273.15)/(2*FARADAY)
	nu = vm/f
	ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}


FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}
