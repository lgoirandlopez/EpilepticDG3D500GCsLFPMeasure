TITLE SK.mod
:SK-type calcium-activated potassium current

: Calcium activated K channel.
: Adpated from Tejada et al. 2014 and initially Adapted from Santhakumar et al. 2005 ( which was a midification from Moczydlowski and Latorre (1983) J. Gen. Physiol. 82

COMMENT
Original Mod Files:
Original names 'gskch.mod' 
Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
https://modeldb.science/155568?tab=2&file=TejadaEtAl2014/bgka.mod
Original modified from Moczydlowski and Latorre (1983) J. Gen. Physiol. 82
Original comment : modified 1/7/2007 by Chris Deister for the GP neuron (to remove some of the background current that existed in Mercer 2007)

Author: Lucas Goirand-Lopez
Modified from Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
NEURON {
	SUFFIX sk
	USEION k READ ek WRITE ik VALENCE 1
	USEION nca READ ncai VALENCE 2
	USEION lca READ lcai VALENCE 2
	USEION tca READ tcai VALENCE 2
	RANGE  gskbar, gsk, ik, qinf, qtau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
	celsius (degC)
	temp = 30 (degC)
	v
	gskbar = 0.001	(mho/cm2)
	cai (mM)
	ncai (mM)
	lcai (mM)
	tcai (mM)
	q10 = 3
}

STATE {	q }

ASSIGNED {
	ik	(mA/cm2)
	gsk	(mho/cm2)
	qinf
	qtau	(ms)
	ek	(mV)
	tadj
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gsk = gskbar*q*q
	ik = gsk*(v-ek)
}
UNITSOFF

INITIAL {
	cai = ncai + lcai + tcai	
	qrate(cai)
	q=qinf
	VERBATIM
	ncai = _ion_ncai;
	lcai = _ion_lcai;
	tcai = _ion_tcai;
	ENDVERBATIM
}

DERIVATIVE state {
	cai = ncai + lcai + tcai
	qrate(cai)
	q' = (qinf - q)/qtau
}

PROCEDURE qrate(cai (mM)) {
	LOCAL alpha, beta
	tadj = q10^((celsius - temp)/10)
	alpha = 1.25e1 * cai * cai
	beta = 0.00025 
:	alpha = 0.00246/exp((12*log10(cai)+28.48)/-4.5)
:	beta = 0.006/exp((12*log10(cai)+60.4)/35)
	
	qinf = alpha/(alpha + beta)
	qtau = 1/(alpha + beta)/tadj
}

UNITSON

