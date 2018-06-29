simprepTemporal <- function(totaltime, regions=NULL, onsets, durations, TR, effectsize, accuracy=0.1, hrf=c("gamma","double-gamma","Balloon"), param=NULL){
  
  if(missing(hrf)){
    hrf <- "double-gamma"
  }
  if(!is.list(onsets)){
    if(is.numeric(onsets)){
      onsets <- list(onsets)
    }else{
      stop("Wrong format of onsets. See ?simprepTemporal for further details.")
    }
  }
  if(!is.list(durations)){
    if(is.numeric(durations)){
      durations <- as.list(durations)
    }else{
      stop("Wrong format of durations. See ?simprepTemporal for further details.")
    }
  }
  if(!is.list(effectsize)){
    if(is.numeric(effectsize)){
      effectsize <- as.list(rep(effectsize, length(onsets)))
    }else{
      stop("Wrong format of effect size. See ?simprepTemporal for further details.")
    }
  }
  if(!is.null(regions)){
    if(regions!=1){
      if((length(onsets)!=regions)||(length(durations)!=regions)){
        stop("Mismatch between number of regions and onsets and/or durations. See ?simprepTemporal for further details.")
      }
    }
  } else {
    regions <- 1
  }
  if(!is.null(param)){
    if(regions!=1){
      if(length(param)!=regions){
        stop("Mismatch between number of regions and param. See ?simprepTemporal for further details.")
      }
    }
  }
  
  temp = list(onsets = onsets,
              durations = durations,
              totaltime = totaltime,
              effectsize = effectsize,
              TR = TR,
              acc = accuracy,
              hrf = hrf,
              param = param)
  
  return(purrr::map(1:regions, function(x) {temp$name <- paste0("region", x); return(temp)}))
}

