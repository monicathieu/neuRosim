data.read.nifti <-
function(filename="", verbose=TRUE){

	require(Rniftilib, quiet=TRUE)
	if(verbose) cat("Reading in 4D nifti file.")
	Y <- nifti.image.read(filename)
	if(verbose){
		cat("done.\n")
		cat("Dimensions: [", dim(Y), "]\n", sep=" ")
	}
	Y <- array(Y, dim=dim(Y))
	return(Y)
}

