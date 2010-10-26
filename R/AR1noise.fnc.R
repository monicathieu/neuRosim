AR1noise.fnc <-
function(dim, sigma, nscan, rho=NULL, anatomical=FALSE, template, verbose=TRUE){

	# dim: vector representing the dimension of the noise image
	# sigma: standard deviation of the noise
	# nscan: number of scans or timepoints
	# rho: autocorrelation parameter
        # anatomical: logical if an anatomical template is used or not
        # template : image array defining the anatomical structure

	if(length(dim)>3){
		stop("Image space with more than three dimensions is not supported.")
	}
	if(is.null(rho)){
		if(verbose==TRUE){
			warning("Default value for autocorrelation parameter is used.")
		}
		rho <- 0.2
	}

	array.dim <- c(dim,rep(1,3-length(dim)),nscan)
        noise.array <- array(0, dim=array.dim)
 
	for(i in 1:nscan){
		if(i==1){
			noise.array[,,,i] <- rnorm(prod(array.dim[-4]), 0, sigma)/sqrt(1-rho^2)
		} else {
			noise.array[,,,i] <- rho * noise.array[,,,i-1] + rnorm(prod(array.dim[-4]), 0, sigma)
		}
	}

	noise <- array(c(noise.array), dim=c(dim, nscan))/sqrt(1-rho^2)

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

