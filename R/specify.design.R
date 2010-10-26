specify.design <-
function(ncond, onsets, durations, T, TR, acc, conv=c("none", "gamma", "double-gamma", "Balloon"), cond.names=NULL, par=NULL){

	# ncond: number of conditions, i.e. columns in the design matrix
	# onsets: vector or matrix representing the onsets in seconds, ncol should match ncond
	# durations: vector or matrix representing the durations in seconds, ncol should match ncond
	# T: duration of the experiment in seconds
	# TR: repetition time in seconds
	# acc: microtime resolution in seconds
	# conv: should the design matrix be convoluted, default is none
	# cond.names: optional names for the conditions
	# par: parameters of the hrf

	if(is.matrix(onsets)){
		if(ncol(onsets)!=ncond){
			stop("Mismatch between number of conditions and onsets.")
		}
	} else if(is.vector(onsets)){
		if(ncond!=1){
			stop("Mismatch between number of conditions and onsets.")
		} else {
			onsets <- matrix(onsets, ncol=ncond)
		}
	}
        if(is.matrix(durations)){
                if(ncol(durations)!=ncond){
                        stop("Mismatch between number of conditions and durations.")
                }
        } else if(is.vector(durations)){
                if(length(durations)==1){
                        durations <- matrix(rep(durations,ncond), ncol=ncond)
                } else if(length(durations)==ncond){
                        durations <- matrix(durations, ncol=ncond)
                } else {
			stop("Mismatch between number of conditions and durations.")
		}
        }
	if(!is.null(par)){
	        if(is.matrix(par)){
        	        if(ncol(par)!=ncond){
                	        stop("Mismatch between number of conditions and par.")
                	}
	        } else if(is.vector(par) && (ncond>1)){
               	        par <- matrix(rep(par,ncond), ncol=ncond)
	       
               	} else {
                        stop("Wrong format of par.")
        	}
	}

	if(missing(conv)){
		conv <- "none"
	}
	if(is.null(cond.names)){
		cond.names <- c(paste("C", 1:ncond, sep=""))
	}

	design.matrix <- matrix(0,(T/TR),ncond)
	ix <- seq(1,T/acc,TR/acc)
	for(i in 1:ncond){
		if(conv=="none"){
	                design.matrix[,i] <- stim.boxcar(T, onsets[,i], durations[,i], acc)[ix]
 		}
		if(conv=="gamma"){
			s <- stim.boxcar(T, onsets[,i], durations[,i], acc)
			s.conv <- convolve(gammaHRF(seq(acc,T,acc), par[,i], verbose=FALSE), rev(s))
			design.matrix[,i] <- s.conv[ix]
		}
		if(conv=="double-gamma"){
                        s <- stim.boxcar(T, onsets[,i], durations[,i], acc)
                        s.conv <- convolve(canonicalHRF(seq(acc,T,acc), par[,i], verbose=FALSE), rev(s))
                        design.matrix[,i] <- s.conv[ix]
 		}
		if(conv=="Balloon"){
			s <- stim.boxcar(T, onsets[,i], durations[,i], acc)
			s.conv <- Balloon(s, T, acc)
			design.matrix[,i] <- s.conv[ix]
		}
	}
	colnames(design.matrix) <- cond.names
	return(design.matrix)
}

