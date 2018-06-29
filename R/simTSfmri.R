#' @importFrom dplyr mutate select starts_with %>%

simTSfmri <- function(design=list(), base=0, nscan=NULL, TR=NULL, SNR=NULL, noise=c("none", "white", "temporal", "low-frequency", "physiological", "task-related", "mixture"), type=c("gaussian","rician"), weights, verbose=TRUE, rho=0.2, freq.low=128, freq.heart=1.17, freq.resp=0.2, vee=1){
  
  if(missing(noise)) noise <- "white"
  
  if(missing(type)) type <- "gaussian"
  
  if(noise=="mixture"){
    if(missing(weights)){
      stop("Weights should be provided with noise=mixture.")
    }
    if(length(weights)!=5){
      stop("Weights vector should have 5 elements.")
    }
    if(sum(weights)!=1){
      stop("The sum of the weights vector should be equal to 1.")
    }
  }
  
  if (is.null(SNR)) stop("SNR not specified, with no default.")
  
  # Might deprecate the case where a design object isn't specified
  if(length(design)==0){
    act <- base
    sigma <- mean(act)/SNR
    if(is.null(TR)){
      stop("TR value is missing.")
    }
    if(is.null(nscan)){
      stop("nscan value is missing.")
    }
  } else if(length(design)>1){
    stop("Multiple regions are undefined for time series.")
  } else {
    convs <- specifydesign(stimfunction(design[[1]]$totaltime, design[[1]]$onsets, design[[1]]$durations, accuracy = design[[1]]$acc),
                           TR = design[[1]]$TR,
                           effectsize = design[[1]]$effectsize,
                           conv = design[[1]]$hrf,
                           param = design[[1]]$param) %>%
      mutate(act = base + rowSums(select(., starts_with("conv."))))
    
    sigma <- mean(convs$act)/SNR
    TR <- design[[1]]$TR
    nscan <- design[[1]]$totaltime/design[[1]]$TR
  }
  
  # Cannot use case_when, regrettably, because it requires all conditions to be evaluable
  # which is not always true
  if (noise == "none") n <- 0
  # why is this wrapped in a c() call? inspect the type of the output of systemnoise()
  # the noise functions output an array with nrow = dim. c() coerces it to a vector when dim = 1
  # TODO: Either remove verbose args from noise functions or add print functionality
  if (noise == "white") n <- c(systemnoise(dim=1, sigma=sigma, nscan=nscan, type=type, verbose=verbose, vee=vee))
  if (noise == "temporal") n <- c(temporalnoise(dim=1, sigma=sigma, nscan=nscan, rho=rho, verbose=verbose))
  if (noise == "low-frequency") n <- c(lowfreqdrift(dim=1, freq=freq.low, nscan=nscan, TR=TR, verbose=verbose))
  if (noise == "physiological") n <- c(physnoise(dim=1, sigma=sigma, nscan=nscan, TR=TR, freq.heart=freq.heart, freq.resp=freq.resp, verbose=verbose))
  if (noise == "task-related") n <- c(tasknoise(act.image=convs$act, sigma=sigma, type=type, vee=vee))
  if (noise=="mixture") {
    # yes, now they're all calculated even ifthe weight was 0; I think this function wasn't that slow to begin with
    # so hopefully the system time issue won't be turrible
    
    n.white <- c(systemnoise(dim=1, sigma=sigma, nscan=nscan, type=type, verbose=verbose, vee=vee))
    n.temp <- c(temporalnoise(dim=1, sigma=sigma, nscan=nscan, rho=rho, verbose=verbose))
    n.low <- c(lowfreqdrift(dim=1, freq=freq.low, nscan=nscan, TR=TR, verbose=verbose))
    n.phys <- c(physnoise(dim=1, sigma=sigma, nscan=nscan, TR=TR, freq.heart=freq.heart, freq.resp=freq.resp, verbose=verbose))
    n.task <- c(tasknoise(act.image=convs$act, sigma=sigma, type=type, vee=vee))
    
    w <- weights
    n <- (w[1]* n.white + w[2]*n.temp + w[3]*n.low + w[4]*n.phys + w[5]*n.task)/sqrt(sum(w^2))
  }
  
  convs <- convs %>%
    mutate(noise = n,
           ts = act + noise - mean(noise)) %>%
    select(-act, -noise)
  
  return(convs)
}

