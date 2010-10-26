spatnoise.fnc <-
function(dim, sigma, nscan, method=c("corr", "gammaRF", "gausRF"), rho=NULL, FWHM=NULL, gamma.shape=NULL, gamma.rate=NULL, anatomical=FALSE, template, verbose=TRUE){

	# dim: vector representing the dimensions of the image space
	# sigma: standard deviation of the noise
	# nscan: number of scans/timepoints
	# method: method defining the spatial correlation structure, default is corr
	# rho: if method=corr, the autocorrelation coefficient
	# FWHM: if method=gammaRF or method=gausRF, the FWHM of the random field kernel
	# gamma.shape: if method=gammaRF, the shape parameter of the gamma distribution
	# gamma.rate: if method=gammaRF, the rate parameter of the gamma distribution
	# anatomical: logical if an anatomical template is used
	# template: array specifying the anatomical structure

        if(length(dim)>3){
                stop("Image space with more than three dimensions is not supported.")
        }
	if(length(dim)==1){
		stop("Spatially correlated noise structures are not defined for vectors.")
	}
	if(missing(method)){
		method <- "corr"
	}
	if(is.null(rho) && method=="corr"){
		if(verbose==TRUE){
			warning("The default value for spatial correlation is used.")
		}
		rho <- 0.75
	}
	if(is.null(FWHM) && (method=="gammaRF" || method=="gausRF")){
		if(verbose==TRUE){
			warning("The default value for the full-width-half-maximum is used.")
		}
		FWHM <- 4
	}
	if(is.null(gamma.shape) && method=="gammaRF"){
		if(verbose==TRUE){
			warning("The default value for the gamma shape parameter is used.")
		}
		gamma.shape <- 6
	}
	if(is.null(gamma.rate) && method=="gammaRF"){
		if(verbose==TRUE){
			warning("The default value for the gamma rate parameter is used.")
		}
		gamma.rate <- 1
	}

	if(method=="corr"){
		if(length(dim)==2){
			noise <- array(0, dim=c(dim,nscan))
			for(z in 1:nscan){
				start <- array(rnorm(prod(dim), 0, sigma), dim=dim)
				noise.scan <- array(0, dim=dim)
				noise.scan[1,1] <- start[1,1]
				for(i in 2:dim[1]){
					noise.scan[i,1] <-rho*noise.scan[(i-1),1] + sqrt(1-rho^2)*start[i,1]
				}
				for(j in 2:dim[2]){
					noise.scan[,j] <- rho*noise.scan[,(j-1)] + sqrt(1-rho^2)*start[,j]
				}
				for(i in 2:dim[1]){
					noise.scan[i,2:dim[2]] <- rho*noise.scan[(i-1),2:dim[2]] + sqrt(1-rho^2)*noise.scan[i,2:dim[2]]
				}
				noise[,,z] <- noise.scan
			}
		} else {
			noise <- array(0, dim=c(dim,nscan))
			for(z in 1:nscan){
				start <- array(rnorm(prod(dim), 0, sigma), dim=dim)
				noise.scan <- array(0, dim=dim)
				noise.scan[1,1,1] <- start[1,1,1]
				for(i in 2:dim[1]){
					noise.scan[i,1,1] <- rho*noise.scan[(i-1),1,1] + sqrt(1-rho^2)*start[i,1,1]
				}
				for(j in 2:dim[2]){
					noise.scan[,j,1] <- rho*noise.scan[,(j-1),1] + sqrt(1-rho^2)*start[,j,1]
				}
				for(k in 2:dim[3]){
					noise.scan[,,k] <- rho*noise.scan[,,(k-1)] + sqrt(1-rho^2)*start[,,k]
				}
				for(i in 2:dim[1]){
					noise.scan[i,2:dim[2],2:dim[3]] <- rho*noise.scan[(i-1),2:dim[2],2:dim[3]] + sqrt(1-rho^2)*noise.scan[i,2:dim[2],2:dim[3]]
				}
				noise[,,,z] <- noise.scan
			}
		}
	        if(anatomical==TRUE){
        	        if(missing(template)){
                	        stop("Template image is required when anatomical==TRUE")
               		}
	                template.time <- array(rep(template,nscan), dim=c(dim,nscan))
       		        ix <- which(template.time!=0)
                	noise[-ix] <- 0
        	}
	}
	if(method=="gausRF"){
		require(AnalyzeFMRI, quietly=TRUE)
		s <- diag(FWHM^2, 3)/(8*log(2))
                if(length(dim)==2){
                        dim.RF <- c(dim,1)
                } else {
                        dim.RF <- dim
                }
		voxdim <- rep(1, length(dim.RF))
		if(anatomical==TRUE){
			if(missing(template)){
				stop("Template image is required when anatomical==TRUE")
			}
			m <- array(ifelse(template!=0), dim=dim.RF)
		} else {
			m <- array(1, dim=dim.RF)
		}
		if(FWHM%%2==0){
			ksize <- FWHM + 1
		} else {
			ksize <- FWHM
		}
	
		noise <- array(0, dim=c(dim.RF,nscan))
		for(z in 1:nscan){
			noise[,,,z] <- array(c(Sim.3D.GRF(d=dim.RF,voxdim=voxdim,sigma=s,ksize=ksize,mask=m,type="field")$mat), dim=dim.RF)
		}
		noise <- array(c(noise), dim=dim)
	}
	if(method=="gammaRF"){
                require(AnalyzeFMRI, quietly=TRUE)
                s <- diag(FWHM^2, 3)/(8*log(2))
                if(length(dim)==2){
                        dim.RF <- c(dim,1)
                } else {
                        dim.RF <- dim
                }
                voxdim <- rep(1, length(dim.RF))
                if(anatomical==TRUE){
                        if(missing(template)){
                                stop("Template image is required when anatomical==TRUE")
                        }
                        m <- array(ifelse(template!=0), dim=dim.RF)
                } else {
                        m <- array(1, dim=dim.RF)
                }
                if(FWHM%%2==0){
                        ksize <- FWHM + 1
                } else {
                        ksize <- FWHM
                }

                noise <- array(0, dim=c(dim.RF,nscan))
                for(z in 1:nscan){
                        n <- Sim.3D.GRF(d=dim.RF,voxdim=voxdim,sigma=s,ksize=ksize,mask=m,type="field")$mat
			gamma.n <- qgamma(pnorm(c(n)), shape=gamma.shape, rate=gamma.rate)
			noise[,,,z] <- array(gamma.n, dim=dim.RF)
                }
                noise <- array(c(noise), dim=dim)
	}

	return(noise)
}

