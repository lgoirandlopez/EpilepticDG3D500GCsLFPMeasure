TITLE lca.mod 
: L-type calcium channel 

COMMENT
Original Mod Files:
Original names 'lca.mod' 
Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
https://modeldb.science/155568?tab=2&file=TejadaEtAl2014/lca.mod
Original modified from Migliore CA3  

Author: Lucas Goirand-Lopez
Adapted from Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX lca
	USEION lca READ elca WRITE ilca VALENCE 2
	USEION ca READ cai, cao VALENCE 2
	RANGE glcabar, einf, etau, glca, elca, ilca
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (milli/liter)

	FARADAY = 96520 (coulomb)
	R = 8.313424  (joule/degC)
}

PARAMETER {
	v (mV)
	celsius (degC)
	glcabar (mho/cm2)
	ki = .001 (mM)
	cai (mM)
	cao (mM)
	tfa=1
}


STATE {
	e
}

ASSIGNED {
	ilca (mA/cm2)
	glca (mho/cm2)
	einf 
	etau (ms)
	elca (mV)  
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	glca = glcabar*e*e*h2(cai)
	ilca = glca * ghk(v, cai, cao)
}

DERIVATIVE states {	: exact when v held constant
	erate(v)
	e' = (einf - e)/etau
}


INITIAL {
	erate(v)
	e = einf
	VERBATIM
	cai=_ion_cai;
	ENDVERBATIM
}

FUNCTION alp(vm(mV)) (1/ms) {
	TABLE FROM -150 TO 150 WITH 200
	alp = 15.69*(-1.0*vm+81.5)/(exp((-1.0*vm+81.5)/10.0)-1.0)
}

FUNCTION bet(vm(mV)) (1/ms) {
	TABLE FROM -150 TO 150 WITH 200
	bet = 0.29*exp(-vm/10.86)
}

PROCEDURE erate(vm(mV)) {
	LOCAL  alpha, beta
	
	alpha = alp(vm)
	beta = bet(vm)
	
	einf = alpha/(alpha+beta)
	etau = 1/(tfa*(alpha+beta))
}


FUNCTION h2(cai (mM)) {
	h2 = ki/(ki+cai)
}

FUNCTION ghk(vm(mV), ci(mM), co(mM)) (mV) {
	LOCAL nu,f

	f = KTF(celsius)/2
	nu = vm/f
	ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {
        KTF = (1000)*R*(celsius +273.15)/(FARADAY) : C/mol-1 E = qV soit J = C.V
}


FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}
