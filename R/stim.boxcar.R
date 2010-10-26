stim.boxcar <-
function(T, onset, duration, acc){

  # T: total time of the design in seconds
  # onset: vector representing the onsets of the stimulus in seconds
  # duration: vector representing the durations of the stimulus in seconds
  # acc: microtime resolution in seconds

  s <- rep(0,T/acc)
  os <- onset/acc
  dur <- duration/acc
  if(length(duration)==1){
	dur <- dur*rep(1,length(onset))
  } else if(length(duration) != length(onset)) {
		stop("Mismatch between number of onset and number of durations.")
  }
  for(i in (1:length(onset))){
	if((os[i]+dur[i]) <= T/acc){
    		s[c(os[i]:(os[i]+dur[i]))] <- 1
	} else {
		s[c(os[i]:(T/acc))] <- 1
	}
  }
  return(s)
}

