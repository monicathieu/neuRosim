data.write.nifti <-
function(V, output.file, check=1, datatype=16){

	require(Rniftilib, quiet=TRUE)
	if(!is.array(V)){
		stop("V is not an array")
	}
	if(missing(output.file)){
		stop("please specify filename for writing nifti image")
	} else if(!is.character(output.file)){
		stop("output.file is not a character")
	}
	IMG <- nifti.image.new()
	IMG$dim <- dim(V)
	check <- nifti.set.filenames(nim=IMG, prefix=output.file, check=check, set_byte_order=1) ##FIXME
	nifti.image.setdatatype(nim=IMG, value=datatype)
	
	if(length(dim(V))==1){
		IMG[] <- V
	} else if(length(dim(V))==2){
		IMG[,] <- V
	} else if(length(dim(V))==3){
		IMG[,,] <- V
	} else if(length(dim(V))==4){
		IMG[,,,] <- V
	} else if(length(dim(V))==5){
		IMG[,,,,] <- V
	} else if(length(dim(V))==6){
		IMG[,,,,,] <-V
	}
	out <- nifti.image.write(IMG)
	invisible(out)
}

