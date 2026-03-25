TITLE BK.mod 
:  BK-type calcium-activated potassium currentl.
: Adpated from Tejada et al. 2014 and initially Adapted from Santhakumar et al. 2005 ( which was a midification from Moczydlowski and Latorre (1983) J. Gen. Physiol. 82

COMMENT
Original Mod Files:
Original names 'CaBK.mod' 
Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
https://modeldb.science/155568?tab=2&file=TejadaEtAl2014/bgka.mod
Original modified from Moczydlowski and Latorre (1983) J. Gen. Physiol. 82

Author: Lucas Goirand-Lopez
Modified from Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)

ENDCOMMENT

NEURON {
	SUFFIX bk
	USEION nca READ ncai VALENCE 2
	USEION lca READ lcai VALENCE 2
	USEION tca READ tcai VALENCE 2
	USEION k READ ek WRITE ik VALENCE 1
	RANGE gbkbar, gbk,  ik, oinf, otau
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(molar) = (1/liter)
	(mM) = (millimolar)		
	
	FARADAY = (faraday) (kilocoulombs)
	R = 8.313424  (joule/degC)
}

CONSTANT {	
	d1 = .84
	d2 = 1.
	k1 = .48e-3	(mM)
	k2 = .13e-6	(mM)
	abar = .28	(/ms)
	bbar = .48	(/ms)
	st=1            (1)
	
}

PARAMETER {
	v (mV)
	celsius (degC)
	gbkbar = .01 (mho/cm2)
	ek (mV)
	cai = 5e-5 (mM)
	lcai		(mV)
	ncai		(mV)
	tcai		(mV)
}

ASSIGNED {
	ik (mA/cm2)
	oinf
	otau (ms)
	gbk (mho/cm2)
}

STATE {
	o   FROM 0 TO 1
}

INITIAL {
	cai = ncai + lcai + tcai
	rate(v,cai)
	o = oinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gbk = gbkbar * o^st
	ik = gbk * (v - ek)
}

DERIVATIVE states {
	rate(v,cai)
	o' = (oinf-o)/otau
}

PROCEDURE rate(vm (mV), c (mM)) {
	LOCAL alpha
	alpha = alp(vm,c)
	otau = 1/(alpha+ bet(vm, c))
	oinf = alpha*otau
}

FUNCTION alp(vm (mV), c (mM)) (1/ms) { :callable from hoc
	alp = c*abar/(c + exp1(k1,d1,vm))
}

FUNCTION bet(vm (mV), c (mM)) (1/ms) { :callable from hoc
	bet = bbar/(1 + c/exp1(k2,d2,vm))
}

FUNCTION exp1(k (mM), d, vm (mV)) (mM) { :callable from hoc
	exp1 = k*exp(-2*d*FARADAY*vm/R/(273.15 + celsius)) 
}

