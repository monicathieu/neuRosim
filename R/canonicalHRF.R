canonicalHRF <-
function(x, par = NULL, verbose=TRUE){
	
	# x: time vector in seconds
	# par: parameter vector of length 5
	#	par[1] : delay of response (relative to onset)
	#	par[2] : delay of undershoot (relative to onset)
	#	par[3] : dispersion of response
	# 	par[4] : dispersion of undershoot
	#	par[5] : scale of undershoot
	# reference HRF function: 
	#	Friston, KJ, Fletcher, P, Josephs, O, Holmes, AP, Rugg, MD and Turner, R (1998). Event-related fMRI: Characterising differential responses. NeuroImage, 7, 30-40.
	# reference default parameter values: 
	#	Glover, GH (1999). Deconvolution of impulse response in event-related BOLD fMRI. NeuroImage, 9, 416-429.

	if(is.null(par)){
		if(verbose==TRUE){
			warning("Default parameters for HRF are used")
		}
        	par <- c(6,12,0.9,0.9,0.35)
	}
        par[6] <- par[1]*par[3]		# time to peak of response
        par[7] <- par[2]*par[4]		# time to peak of undershoot

        (x/par[6])^par[1]*exp(-(x-par[6])/par[3]) - par[5]*(x/par[7])^par[2]*exp(-(x-par[7])/par[4])
}

