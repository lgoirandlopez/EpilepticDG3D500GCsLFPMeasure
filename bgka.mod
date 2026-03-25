TITLE bgka.mod
:Calcium activated K channel.

COMMENT
Original Mod Files:
Original names 'bgka.mod' 
Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
https://modeldb.science/155568?tab=2&file=TejadaEtAl2014/bgka.mod
Original modified from Moczydlowski and Latorre (1983) J. Gen. Physiol. 82

Author: Lucas Goirand-Lopez
Modified from Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX borgka
	USEION k READ ek WRITE ik VALENCE 1
	RANGE gkabar,gka, ik
	GLOBAL ninf,linf,ntau,ltau
	GLOBAL q10, temp, tadj, vmin, vmax
}

PARAMETER {
	gkabar = .01  	(S/cm2)	: 0.12*1e-4 (pS/um2)								

	temp = 30	(degC)		: original temp 
	q10  = 3				: temperature sensitivity

	v 		(mV)
	celsius	(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
	
	vhalfn=-33.6   (mV)
	vhalfl=-83   (mV)
	a0l=0.08      (/ms)
	a0n=0.02    (/ms)
	zetan=-3    (1)
	zetal=4    (1)
	gmn=0.6   (1)
	gml=1   (1)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 		(mA/cm2)
	gka		(S/cm2)
	ek		(mV)
	ninf 		linf
	ntau (ms)	ltau (ms)
	tadj
}
 

STATE { n l }

INITIAL { 
	rates(v)
	n=ninf
	l=linf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gka = gkabar*n*l
	ik = gka*(v-ek)
} 


DERIVATIVE states {    
	rates(v)
	n' = (ninf - n)/ntau
	l' = (linf - l)/ltau

}


PROCEDURE rates(vm) {  

	:TABLE  mtau, minf DEPEND celsius FROM vmin TO vmax WITH 199
	: "m" sodium activation system - act and inact cross at -40
	
	LOCAL a
	tadj = q10^((celsius - temp)/10)
	a = alpn(vm)
	ninf = 1/(1 + a)
	ntau = betn(vm)/(tadj*a0n*(1+a))
	a = alpl(vm)
	linf = 1/(1+ a)
	ltau = betl(vm)/(tadj*a0l*(1 + a))
}


FUNCTION alpn(vm(mV)) {
  alpn = exp(1.e-3*zetan*(vm-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(vm(mV)) {
  betn = exp(1.e-3*zetan*gmn*(vm-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpl(vm(mV)) {
  alpl = exp(1.e-3*zetal*(vm-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betl(vm(mV)) {
  betl = exp(1.e-3*zetal*gml*(vm-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}




