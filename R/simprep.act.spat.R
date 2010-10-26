simprep.act.spat <-
function(nregio, coord, ext=NULL, form=c("cube", "sphere","manual"), fading=FALSE){

        # nregio: number of activated regions
        # coord: matrix of coordinates, the columns represent the xyz-coordinates
        # ext: if form=cube or sphere, the distance between the center and the edge, if form=manual, the number of voxels in each region
        # form: the form of the activated regions
	# fading: logical indicating if the center of the region should be more activated than the edges
	
	if(missing(form)){
		form <- "cube"
	}
	if(!is.matrix(coord)){
		stop("Coord should be a matrix. See ?simprep.act.spat for further details.")
	}
	if(form!="manual"){
		if(nrow(coord)!=nregio){
			stop("Mismatch between rows coord and nregio.  See ?simprep.act.spat for further details.")
		}
	}
	if(is.null(ext)){
		stop("Argument ext is not specified. See ?simprep.act.spat for further details.")
	}
	if(length(ext)==1){
		ext <- rep(ext, nregio)
	} else if(length(ext)!=nregio){
		stop("Mismatch between ext and nregio. See ?simprep.act.spat for further details.")
	}
	if(form=="manual"){
		if(nrow(coord)!=sum(ext)){
			stop("Mismatch between rows coord and ext. See ?simprep.act.spat for further details.")
		}
	}
	if(length(form)==1){
		form <- rep(form, nregio)
	}
		if((ncol(coord)<2) || (ncol(coord)>3)){
		stop("The columns of coord are out of bounds. See ?simprep.act.spat for further details.")
	}
	if(length(fading)==1){
		fading <- rep(fading, nregio)
	}
	
	spat <- list()
	for(i in 1:nregio){
		spat[[i]] <- list()
		name <- paste("region",i, sep="")
		spat[[i]]$name <- name
		spat[[i]]$coord <- coord[i,]
		spat[[i]]$ext <- ext[i]
		spat[[i]]$form <- form[i]
		spat[[i]]$fading <- fading[i]
	}
	return(spat)
}

