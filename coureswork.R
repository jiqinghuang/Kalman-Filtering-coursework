rm(list=ls())

my.kf <- function(v,params){
	
	p <- dim(params$S)[1]
	n <- dim(v)[2]
	eye <- diag(rep(1,p))
	
	S <- params$S
	SV <- params$SV
	SH <- params$SH
	B <- params$B
	A <- params$A
	
	
	h_tt <- matrix(rep(NA,p*n),p,n)
	h_tmt <- matrix(rep(NA,p*n),p,n)
	
	S_tt <- list()
	S_tmt <- list()
	
	# initialisation
	v_t <- v[,1]
	h_tmt[,1] <- params$mu
	S_tmt[[1]] <- S
	K_t <- S_tmt[[1]] %*% t(B) %*%  solve(B %*% S_tmt[[1]] %*% t(B)+SV)
	h_tt[,1] <- (eye - K_t %*% B) %*% h_tmt[,1] + K_t %*% v_t
	S_tt[[1]] <- (eye - K_t %*% B ) %*% S_tmt[[1]]
	
	for (i in 2:n){

	  ###### FILTERING GOES HERE
	  v_t <- v[,i]
	  h_tmt[,i] <- A %*% h_tt[,i-1]
	  S_tmt[[i]] <- A %*% S_tt[[i-1]] %*% t(A) + SH
	  K_t <- S_tmt[[i]] %*% t(B) %*%  solve(B %*% S_tmt[[i]] %*% t(B)+SV)
	  h_tt[,i] <- (eye - K_t %*% B) %*% h_tmt[,i] + K_t %*% v_t
	  S_tt[[i]] <- (eye - K_t %*% B ) %*% S_tmt[[i]]
		
	}
	
	   h_ttau <- matrix(rep(NA,p*n),p,n)
	   S_ttau <- list()
	   S_cov_tau <- list()
	###### SMOOTHING INITIALISATIONS GO HERE
	   h_ttau[ , n] <- h_tt[,n]
  	 S_ttau[[n]] <- S_tt[[n]]
	
  for (j in n:2){
		
		####### SMOOTHING GOES HERE
    J <- S_tt[[j-1]] %*% t(A) %*% solve(S_tmt[[j]])
    h_ttau[ , j-1] <- h_tt[,j-1] + J %*% ( h_ttau[ , j] - h_tmt[, j])
    S_ttau[[j-1]] <- S_tt[[j-1]] + J %*% ( S_ttau[[j]] - S_tmt[[j]]) %*% t(J)
    S_cov_tau[[j-1]] <- S_ttau[[j]] %*% t(J)
	}
	
	h_tt <- h_tt
	h_ttau <- h_ttau
	S_ttau <- S_ttau
	S_cov_tau <- S_cov_tau
	v_hat <- B %*% h_ttau
	return(list(h_tt=h_tt, S_ttau=S_ttau, h_ttau=h_ttau, v_hat=v_hat, S_cov_tau=S_cov_tau))	
}

my.em <- function(v,params,maxiters){
	
  p <- dim(params$S)[1]
  n <- dim(v)[2]
  eye <- diag(rep(1,p))
  
  S <- params$S
  SV <- params$SV
  SH <- params$SH
  B <- params$B
  A <- params$A
  
  
  for (iter in 1:maxiters){
		## update parameters by runnning one EM step
		## will require both filtering and smoothing estimates
		## so will need to call my.kf
	  ## E STEP
    params$SH <- params$SH
    out.kf <- my.kf(v,params)
    h_ttau <- out.kf$h_ttau
    S_ttau <- out.kf$S_ttau
    S_cov_tau <- out.kf$S_cov_tau
    ## M step
    ### update SV once
    params$SV <- 0
    for (i in 1:n){
       update <- B %*% S_ttau[[i]] %*% t(B) + (B %*% h_ttau[,i]-v[,i]) %*% t(B %*% h_ttau[,i]-v[,i])
       params$SV <- params$SV +1/n * update
     }
    ### update SH once
    params$SH <- matrix(0, nrow=2, ncol=2)
    for (i in 2:n){
    first.part <- S_ttau[[i]] - A %*% t(S_cov_tau[[i-1]]) -  S_cov_tau[[i-1]] %*% t(A) + A %*% S_ttau[[i-1]] %*% t(A)
    second.part <- ( h_ttau[ ,i] - A %*% h_ttau[ ,i-1] ) %*%  t( h_ttau[ ,i] - A %*% h_ttau[ ,i-1] )
    params$SH <- params$SH + 1/(n-1) * (first.part + second.part)
    }
    params$SV <- params$SV
    params$SH <- params$SH
  }
	
	
	return(params)
	
}

########## SCRIPT STARTS HERE

v <- t(as.matrix(read.csv('data.csv')))
n <- dim(v)[2]
p <- 2


SVtrue <- 100
SHtrue <- 10*matrix(c(1,0.7,0.7,2),2,2)
theta <- 0.1
params <- list()
params$A <- 0.95*matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
params$B <- matrix(c(1,0),1,2)
params$SV <- SVtrue		# change to identity for part 3
params$SH <- SHtrue		# change to identity for part 3
params$mu <- matrix(c(1,-1),2,1)
params$S <-  diag(rep(1,p))

### PART A and B of QUESTION 2

out.kf <- my.kf(v,params)

par(mfrow=c(3,1))
plot(1:n,v,main='Observed')
lines(1:n,out.kf$v_hat,col='red')
plot(1:n,v,main='Observed')
lines(1:n,out.kf$h_ttau[1,],main='Hidden 1',type='l',col='red')
plot(1:n,v,main='Observed')
lines(1:n,out.kf$h_ttau[2,],main='Hidden 2',type='l',col='red')

## PART C

# initialise parameters
maxiters <- 20
params$SV <- 1
params$SH <- diag(rep(1,p))

params <- my.em(v,params,maxiters)


