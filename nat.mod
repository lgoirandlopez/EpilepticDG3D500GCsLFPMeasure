TITLE nat.mod
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
	SUFFIX nat
	USEION na READ ena WRITE ina VALENCE 1
	RANGE m, h, gnat, gnatbar,vshiftm,vshifth
	RANGE minf, hinf, mtau, htau
	RANGE vshift, vcoef
	GLOBAL q10, temp, tadj, vmin, vmax
}

PARAMETER {
	gnatbar = .12   	(S/cm2)	: 0.12*1e-4 (pS/um2)
	vshiftm = 0	(mV)		: activation voltage shift
	vshifth = 0  (mV)		: inactivation voltage shift 
								

	temp = 6.3	(degC)		: original temp 
	q10  = 3				: temperature sensitivity

	v 		(mV)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
	
	vshift=0 (mV)
	vcoef=1 (mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	gnat		(S/cm2)
	ena		(mV)
	minf 		hinf
	mtau (ms)	htau (ms)
	tadj
}
 

STATE { m h }

INITIAL { 
	mrates(v+vshiftm)
	hrates(v+vshifth)
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gnat = gnatbar*m*m*m*h
	ina = gnat * (v - ena)
} 


DERIVATIVE states {    
	mrates(v+vshiftm)
	hrates(v+vshifth)
	m' = (minf-m)/mtau
	h' = (hinf-h)/htau

}



PROCEDURE mrates(vm (mV)) {  

	:TABLE  mtau, minf DEPEND celsius FROM vmin TO vmax WITH 199
	: "m" sodium activation system - act and inact cross at -40
	
	LOCAL  alpha, beta

	alpha = -0.3*vtrap((vm+60-17+vshift),-5*vcoef)
	beta = 0.3*vtrap((vm+60-45+vshift),5*vcoef)
	
	tadj = q10^((celsius - temp)/10)
	mtau = 1/(alpha+beta)/tadj
	minf = alpha/(alpha+beta)
}


PROCEDURE hrates(vm (mV)) {

	:TABLE  htau, hinf  DEPEND celsius  FROM vmin TO vmax WITH 199
	
	LOCAL  alpha, beta
	
    tadj = q10^((celsius - temp)/10)
	alpha = 0.23/exp((vm+60+5+vshift)/20)
	beta = 3.33/(1+exp((vm+60-47.5+vshift)/-10))
	htau = 1/(alpha+beta)/tadj
	hinf = alpha/(alpha+beta) 
}


FUNCTION vtrap(x (mV),y (mV)) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}



