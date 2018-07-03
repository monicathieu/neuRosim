#' Generate design matrix.
#' 
#' Generates a design matrix to be used as a model for the simulated activation.
#' 
#' @export
#' @importFrom dplyr filter full_join group_by %>%
#' @importFrom tidyr gather nest unnest
#' @importFrom purrr map map2 map_dbl reduce
#' @importFrom stats convolve
#' @importFrom stringr str_sub
#' 
#' @param stim.function A tibble output by \code{\link{stimfunction}} specifying onsets.
#' @param TR Repetition time in seconds.
#' @param effectsize List or number representing the effectsize in each condition.
#' @param conv Should the design matrix be convoluted, default is none.
#' @param param Parameters of the haemodynamic response function.
#' See \code{\link{gammaHRF}} and \code{\link{canonicalHRF}} for more details.
#' @return A tibble with columns specifying the design.
#' 
#' @seealso \code{\link{specifyregion}},\code{\link{gammaHRF}},\code{\link{canonicalHRF}},\code{\link{balloon}}
#' 
#' @examples
#' os <- list(c(20,60),c(15,35))
#' d <- list(20, 10)
#' total <- 80
#' TR <- 2
#' out <- specifydesign(stimfunction(total, os, d, acc = 0.1), TR,
#'   effectsize=list(2,5), conv="double-gamma")
#' \dontshow{rm(os,d,total,TR,out)}
#' 
#' @keywords low-level

specifydesign <- function(stim.function, TR, effectsize, conv=c("none", "gamma", "double-gamma", "Balloon"), param=NULL){
  # overhaul the input structure?
  
  ncond <- sum(startsWith(names(stim.function), "s."))
  

  if(is.list(effectsize)){
    if(length(effectsize)!=ncond){
      stop("Mismatch between number of conditions and effectsize.")
    }
  } else {
    if(is.numeric(effectsize)){
      effectsize <- list(effectsize)
    }else{
      stop("wrong format of effectsize: should be a list or vector")
    }
  }
  if(!is.null(param)){
    if(is.list(param)){
      if(length(param)!=ncond){
        stop("Mismatch between number of conditions and param.")
      }
    } else {
      stop("param should be a list")
    }
  } else {
    # list-ify the null params
    param = vector("list", ncond)
  }
  
  if(missing(conv)){
    conv <- "none"
  }
  
  totaltime <- max(stim.function$microtime)
  
  prep <- stim.function %>%
    gather(key = "cond.name", value = "s", -microtime) %>%
    group_by(cond.name) %>%
    nest(.key = "byMicro") %>%
    mutate(param = map(param, ~.),
           effectsize = map_dbl(effectsize, ~.))

    if (conv=="none") {
      # if no convolution: ... do nothing
      prep <- prep %>%
        mutate(byTR = map(byMicro, ~filter(., microtime %% TR == 0) %>%
                            mutate(conv = s)))
    } else if (conv=="gamma") {
      if(!is.null(param[[1]])){
        prep <- prep %>%
          mutate(conv = map2(byMicro, param, ~convolve(gammaHRF(.x$microtime, FWHM = unlist(.y), verbose = F), rev(.x$s))))
      } else {
        prep <- prep %>%
          mutate(conv = map(byMicro, ~convolve(gammaHRF(.x$microtime, verbose = F), rev(.x$s))))
      }
    } else if (conv=="double-gamma") {
      prep <- prep %>%
        mutate(conv = map2(byMicro, param, ~convolve(canonicalHRF(.x$microtime, param = .y, verbose = F), rev(.x$s))))
    } else if (conv=="Balloon") {
      # this one doesn't use anything specified in param at all? ok...
      prep <- prep %>%
        mutate(conv = map(byMicro, ~balloon(.x$s, totaltime, accuracy)))
    }
    prep <- prep %>%
      # rescale to be be max of 1, then rescale to effect size
      mutate(conv = map2(conv, effectsize, ~.y*(.x/max(.x))),
             byMicro = map2(byMicro, conv, ~mutate(.x, conv = .y)),
             # downsampling occurs here
             byTR = map(byMicro, ~filter(., round(microtime %% TR, 2) == 0)))
  
  return(reduce(prep$byTR, dplyr::full_join, by = "microtime", suffix = str_sub(prep$cond.name, start = 2L)))
}

