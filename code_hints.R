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
	S_tt <- list()
	
	# initialisation
	v_t <- v[,1]
	h_tmt <- params$mu
	S_tmt <- S
	K_t <- S_tmt %*% t(B) %*%  solve(B %*% S_tmt %*% t(B)+SV)
	h_tt[,1] <- (eye - K_t %*% B) %*% h_tmt + K_t %*% v_t
	S_tt[[1]] <- (eye - K_t %*% B ) %*% S_tmt
	
	for (i in 2:n){

		###### FILTERING GOES HERE
		h_tt[,i] <- (eye - K_t %*% B) %*% h_tmt + K_t %*% v_t
		
	}
	
	###### SMOOTHING INITIALISATIONS GO HERE
	
	for (i in n:2){
		
		####### SMOOTHING GOES HERE
		
	}
	
	h_hat <- h_tt
	v_hat <- B %*% h_tt
	return(list(h_hat=h_hat, v_hat=v_hat))	
}

my.em <- function(v,params,maxiters){
	
	for (iter in 1:maxiters){
		## update parameters by runnning one EM step
		## will require both filtering and smoothing estimates
		## so will need to call my.kf
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

plot(1:n,out.kf$h_hat[1,],main='Hidden 1',type='l',col='red')

plot(1:n,out.kf$h_hat[2,],main='Hidden 2',type='l',col='red')

## PART C

# initialise parameters
maxiters <- 20
params$SV <- 1
params$SH <- diag(rep(1,p))

params <- my.em(v,params,maxiters)


