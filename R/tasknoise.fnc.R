tasknoise.fnc <-
function(dim, act.image, sigma){
	
	# dim: vector representing the dimensions of the image space
	# act.image: array defining where and when activation is present
	# sigma: standard deviation of the task-related noise

	if(missing(act.image)){
		stop("An activation array is required.")
	}

	if(length(dim)==1){
		noise <-rnorm(length(act.image), 0, sigma)
	} else {
		noise <- array(rnorm(prod(dim), 0, sigma), dim=dim)
	}
	ix <- which(zapsmall(act.image)!=0)
	noise[-ix] <- 0

	return(noise)
}

