gammaHRF <-
function(x, FWHM=NULL, verbose=TRUE){
	
	# x: time vector in seconds
	# FWHM: Full Width Half Maximum of the Gamma variate function
	# reference:
	#	Buxton, RB, Uludag, K, Dubowitz, DJ and Liu, TT (2004). Modeling the hemodynamic response to brain activation. NeuroImage, 23, S220-S233.

	if(is.null(FWHM)){
		if(verbose==TRUE){
			warning("Default value for FWHM is used.")
		}
		FWHM <- 4
	}
	th <- 0.242*FWHM
	
	1/(th*factorial(3)) * (x/th)^3 * exp(-x/th)
}

