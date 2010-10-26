simprep.act.temp <-
function(T, ncond, nregio=NULL, onsets, durations, TR, acc, hrf=c("gamma","double-gamma","Balloon"), par=NULL){

	# T: duration of the experiment
	# ncond: number of conditions (i.e. columns in the design matrix)
	# nregio: number of regions
	# onsets: vector, matrix or array representing the onsets of the stimulus in seconds
	# durations: vector, matrix or array representing the durations of the stimulus in seconds
	# TR: repetition time in seconds
	# acc: microtime resolution in seconds
	# hrf: Haemodynamic response function (double-gamma is default)
	# par: vector, matrix or array representing the parameters of the haemodynamic response function

	if(missing(hrf)){
		hrf <- "double-gamma"
	}
	if(!is.null(nregio)){
	    if(nregio!=1){
		if(!is.array(onsets) || !is.array(durations)){
			stop("With multiple regions, onsets and/or durations should be an array. See ?simprep.act.temp for further details.")
		} else if((dim(onsets)[3]!=nregio) || (dim(durations)[3]!=nregio)){
			stop("Mismatch between number of regions and onsets and/or durations. See ?simprep.act.temp for further details.")
		}
	    }
	} else {
		nregio <- 1
	}
	if(!is.null(par)){
		if(is.array(par) && (nregio>1)){
			if(dim(par)[3]!=nregio){
				stop("Mismatch between number of regions and par. See ?simprep.act.temp for further details.")
			}
			if(dim(par)[2]!=ncond){
				stop("Mismatch between number of conditions and par. See ?simprep.act.temp for further details.")
			}
		} else if(is.matrix(par) && (ncond>1)){
			if(ncol(par)!=ncond){
				stop("Mismatch between number of conditions and par. See ?simprep.act.temp for further details.")
			}
		} else if(is.vector(par) && (ncond>1)){
			par <- matrix(rep(par,ncond),ncol=ncond)
		} else {
			stop("Wrong format of par. See ?simprep.act.temp for further details.")
		}
	}
	
	temp <- list()
	if(nregio == 1){
		temp[[1]] <- list()
		name <- "region1"
		temp[[1]]$name <- name
		temp[[1]]$ncond <- ncond
		temp[[1]]$onsets <- onsets
		temp[[1]]$durations <- durations
                temp[[1]]$T <- T
		temp[[1]]$TR <- TR
		temp[[1]]$acc <- acc
		temp[[1]]$hrf <- hrf
		temp[[1]]$par <- par 
	} else {
		for(i in 1:nregio){
			temp[[i]] <- list()
			name <- paste("region", i, sep="")
			temp[[i]]$name <- name
			temp[[i]]$ncond <- ncond
			temp[[i]]$onsets <- matrix(onsets[,,i],ncol=ncond)
			temp[[i]]$durations <- matrix(durations[,,i], ncol=ncond)
                        temp[[i]]$T <-T
			temp[[i]]$TR <-TR
			temp[[i]]$acc  <- acc
			temp[[i]]$hrf <- hrf
			if(!is.null(par)){
				temp[[i]]$par <- matrix(par[,,i],ncol=ncond)	
			} else {
				temp[[i]]$par <- NULL
			}
		}
	}
	return(temp)
}

