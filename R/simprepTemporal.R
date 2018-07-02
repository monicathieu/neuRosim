
#' Prepare temporal structure of the data
#' 
#' Prepare a list defining the necessary parameters to specify the temporal structure
#' of the activation data.
#' 
#' @export
#' @importFrom purrr map
#' 
#' @param totaltime Duration of the experiment.
#' @param regions Number of regions. If not specified, it is assumed that all regions have the same
#' design matrix.
#' @param onsets List or vector representing the onsets of the stimulus in seconds.
#' @param durations List or vector representing the durations of the stimulus in seconds.
#' @param TR Repetition time in seconds.
#' @param effectsize List or number representing the effectsize in each condition.
#' @param accuracy Microtime resolution in seconds.
#' @param hrf Haemodynamic response function (double-gamma is default).
#' @param param Vector, matrix or array representing the parameters of the haemodynamic response function.
#' @return A list with the necessary arguments to be used in \code{\link{simVOLfmri}} or \code{\link{simTSfmri}}.
#' 
#' @seealso \code{\link{simVOLfmri}}, \code{\link{simTSfmri}}, \code{\link{simprepSpatial}},
#' \code{\link{specifyregion}}
#' 
#' @examples
#' ncond <- 2
#' os <- list(c(20,60),c(15,35))
#' d <- list(20, 10)
#' effect <- list(7,10)
#' total <- 80
#' TR <- 2
#' out <- simprepTemporal(total, onsets=os, durations=d, TR=TR, 
#'                      effectsize=effect, hrf="double-gamma")
#' \dontshow{rm(os,d,total,TR,out)}
#' 
#' @keywords high-level

simprepTemporal <- function(totaltime, regions=NULL, onsets, durations, TR, effectsize, accuracy=0.1, hrf=c("gamma","double-gamma","Balloon"), param=NULL){
  
  if(missing(hrf)) hrf <- "double-gamma"
  
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
  
  return(map(1:regions, function(x) {temp$name <- paste0("region", x); return(temp)}))
}

