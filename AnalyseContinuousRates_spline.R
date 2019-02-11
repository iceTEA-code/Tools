###############################################################################################################
######################### Analyse Continuous Rates RATES OF RETREAT OR THINNING ###############################

### R version of the Matlab-based iceTEA tool Analyse Continuous Rates (Bayesian penalised spline analysis).

### Performs penalised spline regression for exposure age data, in a horizontal or vertical transect, to derive 
### a profile of retreat/thinning and corresponding rates through time. Both the normally-distributed exposure 
### ages (2 sigma) and sample position uncertainties are used within a Bayesian framework. It is assumed that 
### exposure ages represent either stability or continuous retreat/thinning through time (i.e. no advance/
### thickening occurs).

### Just Another Gibbs Sampler (JAGS) is used to efficiently perform the analysis. This needs to first be
### downloaded (https://sourceforge.net/projects/mcmc-jags/files/) and installed, before running this analysis.

### Written by Niamh Cahill, University College Dublin, and Richard Selwyn Jones, Durham University.

### The paper describing the tool is:
### Jones, R.S., Small, D., Cahill, N., Bentley, M.J. and Whitehouse, P.L., 2019. iceTEA: Tools for plotting 
### and analysing cosmogenic-nuclide surface-exposure data from former ice margins. Quaternary Geochronology.

### USER SPECIFICATIONS ARE REQUIRED WHERE HIGHLIGHTED WITH CAPITAL LETTERS. RESULTS ARE PLOTTED AT THE END.

###############################################################################################################

## Clear workspace - Start fresh ##
rm (list = ls( )) 


#### SET YOUR WORKING DIRECTORY ####
mydir<-"D:/iceTEA" # EXAMPLE DIRECTORY
setwd(mydir)


#### SET THE DATA FILE TO READ FROM YOUR WORKING DIRECTORY ####
# Data should be in 4 columns: Position (e.g. in km for horizontal transect, m for vertical transect),
# Positon error (e.g. in km or m), mean exposure age (years), and age uncertainty (1 sd, years).

data<-read.csv("Example_transect.csv") # EXAMPLE DATA FILE NAME


###############################################################################################################

## Spline functions ##
getknots = function(x, xl = min(x), xr = max(x), nseg = 30, deg = 3){
  # Construct B-spline basis
  dx = (xr - xl) / nseg
  knots = seq(xl - deg * dx, xr + deg * dx, by = dx)
  return(list(knots=knots,dx=dx,deg=3))
}

## Required Libraries ##
library(rjags)
library(R2jags)
library(coda)

## Sort, standardise and set up the data ##
data.sorted<-data[order(data$Age),]
year<-(-data.sorted$Age-mean(-data.sorted$Age))/sd(-data.sorted$Age)
N<-length(year)
position<-(data.sorted$Position-mean(data.sorted$Position))/sd(data.sorted$Position)
sigma.y=data.sorted$PositionError/sd(data.sorted$Position)
sigma.x=(data.sorted$AgeError)/sd(-data.sorted$Age)


################################################## JAGS MODEL #################################################

SplineModel="
model
{
  for(i in 1:N)
  {
  # observed y centered on true y with sd sigma.y
  y[i] ~ dnorm(mu.y[i],tau.y[i]) 
  # observed x is centered on true x with sd sigma.x
  x[i] ~ dnorm(x0[i],tau.x[i]) 
  
  
  # true y is a spline function of true x
  mu.y[i]<-inprod(B.ik[i,],beta.k) 

  # components of the basis function
  for(j in 1:n.knots){
  J[i,j]<-step(x0[i]-knots[j])
  L[i,j]<-((x0[i]-knots[j])^deg)*J[i,j]
  }# end j loop
  
  # uninformative prior for true x
  x0[i]~dnorm(0,0.01)
  
  # Precision parameters (1/var)
  tau.y[i]<-1/(pow(sigma.y[i],2)+pow(nu.y,2))
  tau.x[i]<-pow(sigma.x[i],-2)
  } # end i loop 
  
  nu.y~dunif(0,1)
  # basis functions of x0 (true x)
  B.ik <- pow(-1,(deg + 1)) * (L %*% t(D))
  
  # Spline parameters (coefficients and smoothness)
  for (k in 2:K) {
  beta.k[k]<-delta.k[k]+beta.k[k-1]
  # Delta controls the smoothness
  # Constraint [T(0,)] imposed on delta so that the relative position must decrease over time
  delta.k[k]~dnorm(0,tau.delta)T(,0)
  }# end k loop
  
  ## Priors for spline parameters ##
  beta.k[1] ~ dnorm(0, 0.01) # we need a prior on the first beta (uninformative)
  tau.delta<-pow(sigma.delta,-2)
  sigma.delta ~ dt(0, 2^-2, 1)T(0,)
  beta.d~dnorm(0,0.01)
  
  # Create predictions using xstar
  for(l in 1:n.grid)
  {
  mu.y.pred[l]<-inprod(Bstar.ik[l,],beta.k)
  for(j in 1:n.knots)
  {
  Jstar[l,j]<-step(xstar[l]-knots[j])
  Lstar[l,j]<-((xstar[l]-knots[j])^deg)*Jstar[l,j]
  }
  }# end l loop
  # basis functions for xstar
  Bstar.ik <- pow(-1,(deg + 1)) * (Lstar %*% t(D))
  
}##End model
  "
  
  ## Get required jags data ##
  
  # Because we have x error we need to create the spline basis functions within JAGS
  # Use function getknots to get the knots for the splines
  res= getknots(year,nseg=20) 
  
  ## JAGS data needed for basis functions ##
  n.knots<-length(res$knots)
  D <- diff(diag(n.knots), diff = res$deg + 1) / (gamma(res$deg + 1) * res$dx ^ res$deg)
  K <- dim(D)[1] 
  
  ## Set up prediction data ##
  x_star<-seq(min(year),max(year),length.out=100)
  n.grid=length(x_star)
  
  mydata <- list(N=N,
                 y=position,
                 x=year,
                 xstar=x_star,
                 sigma.y=sigma.y,
                 sigma.x=sigma.x,
                 n.knots=n.knots,
                 deg=res$deg,
                 D=D,
                 knots=res$knots,
                 n.grid=n.grid,
                 K=K) 
  
  line_inits <- list(list("nu.y" = 0.5, "sigma.delta"= 1),
                     list("nu.y" = 0.5, "sigma.delta"= 1),
                     list("nu.y" = 0.5, "sigma.delta"= 1),
                     list("nu.y" = 0.5, "sigma.delta"= 1))

  
  ## Parameters to look at/save ##
  mypars <- c("beta.k","sigma.delta","mu.y.pred")  
  
  
  ## Run the model ##
  run.mcmc <- jags(data=mydata,
                   inits=line_inits,
                   parameters.to.save=mypars, 
                   model.file=textConnection(SplineModel),
                   # MCMC settings
                   n.chains=4, 
                   n.iter=20000, 
                   n.burnin=4000,
                   n.thin=4)
  

  ## Save the list of parameter names ##
  mcmcmodel.jags<-run.mcmc$BUGSoutput$sims.list;
  
  ## Get the parameter array to check convergence plots/parameters ##
  mcmc.array<-run.mcmc$BUGSoutput$sims.array 

  
  ## Sort outputs ##
  ypred=((mcmcmodel.jags$mu.y.pred)*sd(data.sorted$Position))+mean(data.sorted$Position)
  x_arr=-((x_star*sd(-data.sorted$Age))+mean(-data.sorted$Age))
  
  medy<-apply(ypred,2,median)
  posu68<-apply(ypred,2,quantile,prob=0.84,na.rm=T)
  posl68<-apply(ypred,2,quantile,prob=0.16,na.rm=T)
  posl68bck<-posl68[length(posl68):1]
  posu95<-apply(ypred,2,quantile,prob=0.975,na.rm=T)
  posl95<-apply(ypred,2,quantile,prob=0.025,na.rm=T)
  posl95bck<-posl95[length(posl95):1]
  posxbck<-x_arr[length(posl95):1]
  posxxx<-c(x_arr,posxbck)
  posyyy68<-c(posu68,posl68bck)
  posyyy95<-c(posu95,posl95bck)
  
  
  ## Compute rates ##
  y_diff=apply(ypred,1,diff)
  x_diff = diff(x_arr)
  rate = apply(y_diff, 2, function(x) x/x_diff)
  
  medrate<-apply(rate,1,median)
  ratel68<-apply(rate,1,quantile,probs=0.16)
  rateu68<-apply(rate,1,quantile,probs=0.84)
  ratel95<-apply(rate,1,quantile,probs=0.025)
  rateu95<-apply(rate,1,quantile,probs=0.975)
  x_arr_rate<-x_arr[-1]
  ratexbck<-x_arr_rate[length(ratel95):1]
  ratexxx<-c(x_arr_rate,ratexbck)
  ratel68bck<-ratel68[length(ratel68):1]
  ratel95bck<-ratel95[length(ratel95):1]
  rateyyy68<-c(rateu68,ratel68bck)
  rateyyy95<-c(rateu95,ratel95bck)

  
  
################################################# PLOT RESULTS ################################################

  color_transparent.one <- rgb(0,0,0,alpha=0.1)
  par(cex=1, mai=c(0.3,0.3,0.1,0.1))
  
  ## Plot modelled spline profile ##
  par(fig=c(0.1,0.98,0.4,1))
  plot(1, type="n", ylab='', xlab='', xaxt="n", ylim=range(data.sorted$Position,posyyy95),xlim=rev(range(data.sorted$Age/1000))) #pin=c(3,3)
  grid (NULL,NULL, lty = 'solid', col = "snow2")
  polygon(posxxx/1000, posyyy95, col=color_transparent.one, border = NA)
  polygon(posxxx/1000, posyyy68, col=color_transparent.one, border = NA)
  lines(x_arr/1000,medy,col="black",lwd=2)
  ageupp<-data.sorted$Age+data.sorted$AgeError
  agelow<-data.sorted$Age-data.sorted$AgeError
  arrows(agelow/1000, data.sorted$Position, ageupp/1000, data.sorted$Position,col="rosybrown2",lwd=4,length=0)
  points(data.sorted$Age/1000,data.sorted$Position,pch=21,col="red2",bg="rosybrown2",lwd=2)
  mtext("Relative position", side=2, line=3)

  ## Plot modelled rates
  par(fig=c(0.1,0.98,0.1,0.41),new=TRUE)
  plot(2, type="n",xlab='', ylab='', ylim=range(medrate,rateyyy95),xlim=rev(range(data.sorted$Age/1000)))
  grid (NULL,NULL, lty = 'solid', col = "snow2")
  polygon(ratexxx/1000, rateyyy95, col=color_transparent.one, border = NA)
  polygon(ratexxx/1000, rateyyy68, col=color_transparent.one, border = NA)
  lines(x_arr_rate/1000,medrate,col="black",lwd=2)
  mtext("Exposure age (ka)", side=1, line=3)
  mtext("Rate (position unit / year)", side=2, line=3)

  
  