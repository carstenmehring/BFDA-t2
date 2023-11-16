## ======================================================================
## These are the simulation files for a paired t-test
## ======================================================================

# ---------------------------------------------------------------------
# sample.function: return a data frame/ matrix / vector with simulated raw data

# get two samples with a specific standardized mean differences
# ES = dz (standardized difference scores)
sample.t2.paired <- function(n, ES, options.sample=NULL) {	

	if (options.sample[[1]]=='NE') {
		l1 <- 1
		l2 <- 1
	} else if (options.sample[[1]]=='SE') {
		l1 <- 1
		l2 <- -1
	} else if (options.sample[[1]]=='SW') {
		l1 <- -1
		l2 <- -1
	} else if (options.sample[[1]]=='NW') {
		l1 <- -1
		l2 <- 1
	} else if (options.sample[[1]]=='C') {
		X <- cbind( c(-1,0,1), c(-1,0,1) )
		a <- sample(1:3,size=n,replace=TRUE)
		l1 <- X[a,1]
		l2 <- X[a,2]
	}
	
	x <- rnorm(n, l1*ES, sd=1)
	y <- rnorm(n, l2*ES, sd=1)
	
	return(cbind(x, y))
}

# ---------------------------------------------------------------------
# select.function: select a specified number of accumulating data from 
# the data frame/ matrix that was simulated with sample.function

select.t2.paired <- function(MAXSAMP, n) {
	return(MAXSAMP[1:n, ])
}


# ---------------------------------------------------------------------
# freq.test.function: return p.value, test statistic, and empirical ES

freq.test.t2.paired <- function(SAMP, alternative=NULL, options.sample=NULL) {

	t1 <- t.test(SAMP[,1], mu=0, alternative=alternative)

	# see http://journal.frontiersin.org/article/10.3389/fpsyg.2013.00863/full

	# must returns these values
	return(list(
		statistic = t1$statistic,
		p.value = t1$p.value,
		emp.ES = t1$statistic / sqrt(length(SAMP))
	))
}

# ---------------------------------------------------------------------
# Check definition of prior

prior.check.t2.paired <- function(prior=NULL){
  
  if(!is.list(prior)){
    if(!is.null(prior) == TRUE) {
      stop("Argument prior needs to be specified as a list.")
    } else {
      prior <- list("Cauchy", list(prior.location = 0, prior.scale = sqrt(2)/2))
    }}
  
  match.arg(prior[[1]], c("Cauchy", "t", "normal"))
  
  if(prior[[1]] %in% c("Cauchy", "t")){
    if(is.null(prior[[2]][["prior.location"]])){
      warning("Prior location not defined. Default specification will be used.")
      prior[[2]][["prior.location"]] <- 0
    } 
    if(is.null(prior[[2]][["prior.scale"]])){
      warning("Prior scale not defined. Default specification will be used.")
      prior[[2]][["prior.scale"]] <- sqrt(2)/2
    }
  }
  
  if(prior[[1]] == "Cauchy"){
    if(any(!is.element(names(prior[[2]]), c("prior.location", "prior.scale")))){
      warning("Cauchy distribution only takes parameters prior.location and prior.scale.")
    }
  } else if (prior[[1]] == "t"){
    if(any(!is.element(names(prior[[2]]), c("prior.location", "prior.scale", "prior.df")))){
      warning("t distribution only takes parameters prior.location, prior.scale, and prior.df.")
    }
    if(is.null(prior[[2]][["prior.df"]])) {
      warning("Prior degrees of freedom not defined. Default specifications will be used.")
      prior[[2]][["prior.df"]] <- 1
    }
  } else if (prior[[1]] == "normal") {
    if(any(!is.element(names(prior[[2]]), c("prior.mean", "prior.variance")))){
      warning("Normal distribution only takes parameters prior.mean and prior.variance.")
    }
    if(is.null(prior[[2]][["prior.mean"]])){
      warning("Prior mean not defined. Default specification will be used.")
      prior[[2]][["prior.mean"]] <- 0
    } 
    if(is.null(prior[[2]][["prior.variance"]])){
      prior[[2]][["prior.variance"]] <- 1
      warning("Prior variance not defined. Default specification will be used.")
    }
  }
  return(prior)
}


# ---------------------------------------------------------------------
# BF.test.function: return log(BF10)


BF.test.t2.paired <- function(SAMP, alternative=NULL, freq.test=NULL, prior=NULL, ...) {

	ones <- rep(1,nrow(SAMP))
	d1=SAMP[,1];
	d2=SAMP[,2];
	mlm <- lm(cbind(d1,d2) ~ -1 + ones)
	BFmlm <- BF(mlm,hypothesis="ones_on_d1>0 & ones_on_d2<0")

	# returns the log(BF10)

	BF10<-BFmlm$BFmatrix_confirmatory[1,2]
	if (is.infinite(BF10) || BF10>1e20) {
		BF10<-1e20
	}
	if (BF10<0) {
		return(NaN)
	}
	if (BF10<1e-20) {
		BF10<-1e-20
	}
	
	return(as.numeric(log(BF10)))

}
