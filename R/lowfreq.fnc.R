lowfreq.fnc <-
function(dim, freq=NULL, nscan, TR, anatomical=FALSE, template, verbose=TRUE){
	# dim: vector representing the dimensions of the image space
	# freq: frequency of low-frequency drift in Hz (default=128)
	# nscan: number of scans or timepoints
	# TR: repetition time
	# anatomical: logical if an anatomical template is used or not
	# template: image array defining the anatomical structure
	spm_drift <- function(N, K) {
    	  n <- 0:(N-1)
    	  C <- matrix(0, nrow=N, ncol=K)

     	  C[,1] <- 1/sqrt(N)
    	  for(k in 2:K) {
       	  C[,k] = sqrt(2/N)*10*cos(pi*(2*n+1)*(k-1)/(2*N))
    	  }
	  return(C)
	}


	if(is.null(freq)){
		if(verbose==TRUE){
			warning("Default frequency is used.")
		}
		freq <- 128
	}

	n <- floor( 2*(nscan*TR)/freq + 1 )
	if(n<3){
		stop("Number of basis functions is too low. Lower frequency or longer scanning time should be provided.")
	}
	drift.base <- spm_drift(nscan, n)[,-1]  
        drift.image <- array(rep(1, prod(dim)), dim=c(dim))
  	drift.array <- drift.image %o% rowSums(drift.base)

        if(anatomical==TRUE){
                if(missing(template)){
                        stop("Template image is required when anatomical==TRUE")
                }
                template.time <- array(rep(template,nscan), dim=c(dim,nscan))
                ix <- which(template.time!=0)
                drift.array[-ix] <- 0
        }

	return(drift.array)
}

