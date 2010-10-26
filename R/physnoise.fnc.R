physnoise.fnc <-
function(dim, sigma, nscan, TR, freq.heart=NULL, freq.resp=NULL, anatomical=FALSE, template, verbose=TRUE){

	# dim: vector representing the dimensions of the image space
	# sigma: standard deviation of the noise
	# nscan: number of scans or timepoints
	# freq.heart: frequency of the heart beat related noise
	# freq.resp: frequency of the respiratory rate related noise
        # anatomical: logical if an anatomical template is used or not
        # template : image array defining the anatomical structure

        if(length(dim)>3){
                stop("Image space with more than three dimensions is not supported.")
        }
	if(is.null(freq.heart)){
		if(verbose==TRUE){
			warning("Default value for heart beat frequency is used.")
		}
		freq.heart <- 1.17
	}
	if(is.null(freq.resp)){
		if(verbose==TRUE){
			warning("Default value for respiratory rate frequency is used.")
		}
		freq.resp <- 0.2
	}

	HB <- 2*pi*freq.heart*TR
	RR <- 2*pi*freq.resp*TR
	t <- 1:nscan

	HRdrift <- sin(HB*t) + cos(RR*t)
	sigma.HR <- sd(HRdrift)
	HRweight <- sigma/sigma.HR

	noise <- array(rnorm(prod(dim)*nscan, 0, 1), dim=c(dim, nscan)) + HRweight*HRdrift

        if(anatomical==TRUE){
                if(missing(template)){
                        stop("Template image is required when anatomical==TRUE")
                }
                template.time <- array(rep(template,nscan), dim=c(dim,nscan))
                ix <- which(template.time!=0)
                noise[-ix] <- 0
        }

       return(noise)
}

