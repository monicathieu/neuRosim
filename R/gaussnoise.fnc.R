gaussnoise.fnc <-
function(dim, sigma, nscan, anatomical=FALSE, template, verbose=TRUE){
	
	# dim: vector representing the dimensions of the noise image
	# sigma: standard deviation of the noise
	# nscan: number of scans or timepoints
	# anatomical: logical if an anatomical template is used or not
	# template : image array defining the anatomical structure

	if(length(dim)>3){
		stop("Image space with more than three dimensions is not supported.")
	}

	noise <- array(rnorm(prod(dim)*nscan, 0, sigma), dim=c(dim,nscan))

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

