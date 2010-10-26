sim.restingstate.ts <- function(nscan, TR, SNR, noise=c("none", "white", "temporal", "low-frequency", "physiological", "task-related", "mixture"), temp=c("AR(1)", "state-space"), weights, verbose=TRUE, rho=NULL, lowfreq=NULL, heartfreq=NULL, respfreq=NULL){
# Synthesizes a single time series x representing resting state activity.
# The fluctuation frequencies f are limited to a square passband
# 0.01 Hz <= f <= 0.1 Hz. TR is the repetition time (needed to compute
# the passband limits), expressed in seconds. N is the required number of
# samples (needs not be a power of 2).
#
# References:
# [1] C.G. Fox, Computers & Geoscience, Vol. 13, pp. 369-374, 1987.
# [2] M. Fukunaga, Magnetic Resonance Imaging, Vol. 24, pp. 979-992, 2006.
	if(missing(noise)){
		noise <- "white"
	}
	if(noise=="temporal" || noise=="mixture"){
		if(missing(temp)){
			temp <- "AR(1)"
		}
	}
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

fl <- 0.01                         # Lower limit of the resting state passband (Hz).
fu <- 0.1                          # Upper limit
Ampl <- 2                          # Maximum percentage BOLD signal change
N <- nscan
f <- 1/TR                          # Sampling frequency
NOS <- 2^ceiling(log2(abs(N)))     # Effective number of samples is the next power of 2 following N.
Nl <- round(fl*NOS/f)              # Express frequencies as fractions of the number of samples.
Nu <- round(fu*NOS/f) 

## Spectrum construction
Z <- rep(0,length=NOS) 
Phase <- 2*pi*runif(NOS/2-1,min=0, max=1)-pi                            # Randomize the phases. 
z <- rep(0,length=NOS/2-1) 
z[Nl:Nu] <- complex(real=cos(Phase[Nl:Nu]),imaginary=sin(Phase[Nl:Nu])) # Construct the passband and force Hermitian.
Z[2:(NOS/2)] <- z                                                       # Lower spectrum half
Z[(NOS/2+2):NOS] <- Conj(z[length(z):1])                                # Upper half

x <- Re(fft(Z,inverse=TRUE)/NOS)   # Inverse Fourier transform (output should be real).
x <- Ampl*x/max(x)                 # Normalize x.
x <- x[1:N]                        # Retain only the first N samples (slightly distorts the power spectrum).

  sigma <- sd(x)/SNR

	if(noise=="none"){
		n <- 0
	}
	if(noise=="white"){
		n <- c(gaussnoise.fnc(c(1), sigma, nscan, verbose=verbose))
	}
	if(noise=="temporal"){
		if(temp=="AR(1)"){
			n <- c(AR1noise.fnc(c(1), sigma, nscan, rho, verbose=verbose))
		}
		if(temp=="state-space"){
			stop("State-space temporal noise model is not yet implemented.")
		}
	}
	if(noise=="low-frequency"){
		n <- c(lowfreq.fnc(c(1), lowfreq, nscan=nscan, TR=TR, verbose=verbose))
	}
	if(noise=="physiological"){
		n <- c(physnoise.fnc(c(1), sigma, nscan, TR, heartfreq, respfreq, verbose=verbose))
	}
	if(noise=="task-related"){
		stop("Task noise is not defined for resting state data.")
	}
	if(noise=="mixture"){
		if(weights[1]==0){
			n.white <- 0
		} else {
			n.white <- c(gaussnoise.fnc(c(1), sigma, nscan, verbose=verbose))
		}
		if(weights[2]==0){
			n.temp <- 0
		} else {
                	if(temp=="AR(1)"){
                        	n.temp <- c(AR1noise.fnc(c(1), sigma, nscan, rho, verbose=verbose))
                	}
                	if(temp=="state-space"){
                        	stop("State-space temporal noise model is not yet implemented.")
                	}
		}
		if(weights[3]==0){
			n.low <- 0
		} else {
			n.low <- c(lowfreq.fnc(c(1), lowfreq, nscan=nscan, TR=TR, verbose=verbose))
		}
		if(weights[4]==0){
			n.phys <- 0
		} else {
			n.phys <- c(physnoise.fnc(c(1), sigma, nscan, TR, heartfreq, respfreq, verbose=verbose))
		}
		if(weights[5]==0){
			n.task <- 0
		} else {
			stop("Task noise is not defined for resting state data.")
		}
		w <- weights
		n <- w[1]* n.white + w[2]*n.temp + w[3]*n.low + w[4]*n.phys + w[5]*n.task
	}
	ts <- x + n
	return(ts)
}
