TITLE nap.mod 
: initially persistent sodium current for nucleus accumbens 

COMMENT
Author: Lucas Goirand-Lopez adapted from
Lee J. (2007). Fast Rhythmic Bursting Cells: The Horizontal Fiber System in the Cat’s Primary Visual Cortex Penn McNair Research Journal. 1(1)

Original Mod Files:
Original names 'nap.mod' 
Lee J. (2007). Fast Rhythmic Bursting Cells: The Horizontal Fiber System in the Cat’s Primary Visual Cortex Penn McNair Research Journal. 1(1)
https://modeldb.science/125857?tab=2&file=FRB/nap.mod
Mod files modified from Moczydlowski and Latorre (1983) J. Gen. Physiol. 82

Original comments :

Magistretti J, Alonso A (1999). "Biophysical properties and slow
voltage-dependent inactivation of a sustained sodium current in entorhinal
cortex layer-II principal neurons." J Gen Phys, 114: 491-509.

Traub RD, Buhl EH et al (2003). "Fast rhythmic bursting can be induced in
layer 2/3 cortical neurons by enhancing persistent na+ conductance or by
blocking BK channels." J Neurophys 89: 909-921.

Jason Moyer 2004 - jtmoyer@seas.upenn.edu
ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S)  = (siemens)
}
 
NEURON {
	SUFFIX nap2
	USEION na READ ena WRITE ina
	RANGE m, gnap, gnapbar, mvhalf, mslope, hvhalf, hslope
	RANGE taumcoef, minf, hinf
		
	GLOBAL q10
}
 
PARAMETER {
	gnapbar   =   4e-5 (S/cm2)	: 4e-5 in soma; 1.3802e-7 in dends

	mvhalf = -52.6		(mV)	: Magistretti 1999, Fig 4
	mslope = -4.6		(mV)	: Magistretti 1999, Fig 4

	hvhalf = -48.8		(mV)	: Magistretti 1999, Fig 4
	hslope = 10.0		(mV)	: Magistretti 1999, Fig 4
	
	temp = 6.3	(degC)		: original temp 
	q10 = 3
	taumcoef = 1
}
 
STATE { m h }
 
ASSIGNED {
	ena		(mV)
	v 		(mV)
	ina		(mA/cm2)
	gnap		(S/cm2)

	minf
	hinf	

	taum	(ms)			: Traub 2003, Table A2
	tadj
   }
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gnap = gnapbar * m * h  
        ina = gnap * ( v - ena )
:        VERBATIM
:        	printf("Ena is %g\n", ena);
:        ENDVERBATIM
}
 

 
INITIAL {
	rates(v,mvhalf,mslope,hvhalf,hslope,taumcoef)
	
	m = minf
	h = hinf
}

:FUNCTION_TABLE tauh(v(mV))  (ms)		: Magistretti 1999, Fig 8A

DERIVATIVE state { 
        rates(v,mvhalf,mslope,hvhalf,hslope,taumcoef)
        m' = (minf - m) / (taum/tadj)
        h' = (hinf - h) / (tauh(v)/tadj)    
}
 
PROCEDURE rates(v (mV), mvhalf (mV), mslope (mV), hvhalf (mV), hslope (mV), taumcoef) {  

	minf = 1 / (1 + exp( (v - mvhalf) / mslope))
	hinf = 1 / (1 + exp( (v - hvhalf) / hslope))
	
	tadj = q10^((celsius - temp)/10)
	
	UNITSOFF
	if (v < -40) {			: Traub 2003, Table A2
		taum = taumcoef*(0.025 + 0.14 * exp( (v + 40 ) / 10))
	} else {
		taum = taumcoef*(0.02 + 0.145 * exp( (-v - 40) / 10))
	}
		
	UNITSON
}

FUNCTION tauh(vm (mV)) {  :Traps for 0 in denominator of rate eqns.
        tauh = 3920*exp(-.5*((vm+68)/29)^2) + 2000
}
 
 
