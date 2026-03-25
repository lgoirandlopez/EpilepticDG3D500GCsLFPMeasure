TITLE gabaal.mod 
:Tonic GABAA leak, adapted from Yim et al (2015))

COMMENT

Original Mod File:
Original name 'ichan2.mod'
Santhakumar V, Aradi I, Soltesz I (2005) J Neurophysiol 93:437-53 
https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=51781&file=%2fdentategyrusnet2005%2fichan2.mod
Morgan RJ, Soltesz I (2008) Proc Natl Acad Sci U S A 105:6179-84
Morgan RJ, Santhakumar V, Soltesz I (2007) Prog Brain Res 163:639-58
Cutsuridis V, Cobb S, Graham BP (2009) Hippocampus 20(3):423-46 

Current version by A. Hanuschkin <AH, 2011> for:
Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

 - A tonic (leak) GABAA conductance to be modified during epilepsy (see Young CC, Stegen M, Bernard R, Muller M, Bischofberger J, Veh RW, Haas CA, Wolfart J (2009) J Physiol 587:4213-4233)

Mod File history:

I_GABAA: (tonic GABAA leak (see above), added in Yim et al (2015))
* replicated from I_leak
 
ENDCOMMENT

:INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX igabaal
	NONSPECIFIC_CURRENT igabaa
	RANGE  ggabaa, igabaa, egabaa
	GLOBAL q10, temp, vmin, vmax
}

PARAMETER {
	ggabaa =  0.722e-05	(mho/cm2)    				: GABAA (gbar and reversal poti)
 	egabaa = -70	(mV)
	
	: Temperature dependence
	temp = 6.3	(degC)		: original temp 
	q10  = 3				: temperature sensitivity

	celsius	= 6.3	(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(um) = (micron)
} 

ASSIGNED {		
	v (mV) 
	
	igabaa (mA/cm2)					: GABAA 
} 

BREAKPOINT {
	igabaa = ggabaa*(v-egabaa)
}

INITIAL {
	VERBATIM
	return 0;
	ENDVERBATIM
}