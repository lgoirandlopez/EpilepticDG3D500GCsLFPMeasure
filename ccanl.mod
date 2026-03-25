TITLE ccanl.mod
:Calcium concentration given inward calcium current and calcium removal

COMMENT 
Original Mod Files:
Original names 'ccanl.mod' 
Tejada, J., Garcia-Cairasco, N., & Roque, A. C. (2014). PLoS Computational Biology, 10(5)
https://modeldb.science/155568?tab=2&file=TejadaEtAl2014/bgka.mod

Old text :
calcium accumulation into a volume of area*depth next to the
membrane with a decay (time constant tau) to resting level
given by the global calcium variable cai0_ca_ion

Warning by Ted Carnevale 2015 (in Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308 model):
The expression that this mechanism uses to calculate the contribution of ica to the rate of change of calcium concentration in the shell is 
-ica*(1e7)/(depth*FARADAY)
but it should really be
-ica*(1e7)/(depth*2*FARADAY)
because the valence of ca is 2.  The result of this omission is that the mechanism behaves as if the shell is only 1/2 as thick as the value specified by the depth parameter.
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}

NEURON {
	SUFFIX ccanl
	USEION nca READ ncai, inca, enca WRITE enca, ncai VALENCE 2
	USEION lca READ lcai, ilca, elca WRITE elca, lcai VALENCE 2
	USEION tca READ tcai, itca, etca WRITE etca, tcai VALENCE 2
	RANGE caiinf, catau, cai, ncai, lcai,tcai, eca, elca, enca, etca
}

UNITS {
	(mV) = (millivolt)
	(molar) = (1/liter)			: moles do not appear in units
	(mM)	= (milli/liter)
	(mA)	= (milliamp)
	FARADAY = 96520 (coulomb)
	R = 8.313424 (joule/degC)
}


PARAMETER {
	celsius = 6.3 (degC)
	depth = 200	(nm)		: depth of shell
	catau	= 10	(ms)		: rate of calcium removal, changed from 200 to 80 (H.Markram)
	caiinf	= 50e-6 (mM)
	cao = 2 (mM)
	ica (mA/cm2)
	inca (mA/cm2)
	ilca (mA/cm2)
	itca (mA/cm2)
	cai= 50.e-6 (mM)
}

ASSIGNED {
	enca (mV)
	elca (mV)
	etca (mV)
	eca (mV)
}

STATE {
	ncai (mM)
	lcai (mM)
	tcai (mM)
}

INITIAL {
	VERBATIM
	ncai = _ion_ncai;
	lcai = _ion_lcai;
	tcai = _ion_tcai; 
	ENDVERBATIM
	ncai=caiinf/3
	lcai=caiinf/3
	tcai=caiinf/3
	cai = caiinf
	eca = ktf() * log(cao/caiinf)
	enca = eca
	elca = eca
	etca = eca
}
	
BREAKPOINT {
	SOLVE integrate METHOD derivimplicit
	cai = ncai+lcai+tcai	
	eca = ktf() * log(cao/cai)	
	enca = eca
	elca = eca
	etca = eca
}

DERIVATIVE integrate { 
ncai' = -(inca)/depth/FARADAY * (1e7) + (caiinf/3 - ncai)/catau
lcai' = -(ilca)/depth/FARADAY * (1e7) + (caiinf/3 - lcai)/catau
tcai' = -(itca)/depth/FARADAY * (1e7) + (caiinf/3 - tcai)/catau
}

FUNCTION ktf() (mV) {
	ktf = (1000)*R*(celsius +273.15)/(2*FARADAY)
} 
