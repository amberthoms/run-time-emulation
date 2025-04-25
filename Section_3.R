
# Section 3: Designing Run Locations --------------------------------------

# Figures 10-16

# Functions ---------------------------------------------------------------
##############################################################################################.
### Define simple Bayes Linear emulator for single input ###
simple_BL_emulator_v1 <- function(x,              # the emulator prediction point
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths     (1 if not specified otherwise)
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])     (1 if not specified otherwise)
                                  E_f = 0         # prior expectation of f: E(f(x)) = 0     (0 if not specified otherwise)
){
  
  # store length of runs D  
  n <- length(D) #different for different no. runs of expensive model so must define n like this
  
  ### Define Covariance structure of f(x): Cov[f(x),f(xdash)] ###
  Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-(x-xdash)^2/theta^2) #specific to our f
  
  
  ### Define 5 objects needed for BL adjustment ### #E[f(x)],E[D],Var[f(x)],Var[D],Cov[f(x),D]
  # Create E[D] vector
  E_D <- rep(E_f,n) #repeat the prior expectation n times
  
  # Create Var_D matrix:
  Var_D <- matrix(0,nrow=n,ncol=n) #nxn zero matrix
  for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i],xD[j]) #Cov fn of x_D elements (as Cov to do with distance between x values)
  
  # Create E[f(x)]
  E_fx <- E_f #prior given in input
  
  # Create Var_f(x) 
  Var_fx <- sigma^2 #sigma ONLY given in input so MUST SQUARE
  
  # Create Cov_fx_D row vector
  Cov_fx_D <- matrix(0,nrow=1,ncol=n) #1xn Cov[f(x),D] #rows=#rows f(x)=1, #columns=#rows D=n
  for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j]) #Cov fn has x_D elements as inputs
  
  
  ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
  ED_fx   <-  E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)   
  VarD_fx <-  Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)  
  
  ### return emulator expectation and variance ###
  return(c("ExpD_f(x)"=ED_fx,"VarD_f(x)"=VarD_fx))  #"..." are the column names!
  
}
### End simple Bayes Linear emulator for single input ###
##############################################################################################.




# Figure 10 ---------------------------------------------------------------
## grids bad because can miss periodicity ----------------------------------
badgrid_plot_BL_emulator_V1 <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))^T
    maintitle=NULL  # title of the plot, blank if NULL
){
  
  ### plot emulator output ###
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-2,2),ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=maintitle) # main argument: plot title (INPUT TO SPECIFY)
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  
  ### plot true function: we would not normally be able to do this! ###
  lines(xP,f(xP),lwd=2,lty=1)
  
  ### Plot the runs ### 
  points(xD,D,pch=21,col=1,bg="green",cex=1.5) #PLOT THE TRUE RUNS OF THE EXPENSIVE MODEL
  legend('bottomright',legend=c("Emulator Expectation","Emulator Prediction Interval",
                                "True function f(x)","Model Evaluations"),
         lty=c(1,1,1,NA),pch=c(NA,NA,NA,16),col=c("blue","red",1,"green"),lwd=2.5,pt.cex=1.3)
}



f <- function(x) cos(10*pi*x)

xD <- seq(0,1,0.2) #evenly spaced

D <- f(xD)
xP <- seq(0.001,0.999,len=201)
#em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))
#must make theta smaller

#png(file="grid_misses_periodicity.png",width=2400,height=1200,res=200)

em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))
badgrid_plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Emulator Output: Limitations of Grid Designs")

#dev.off()





# Figure 11 ---------------------------------------------------------------

all_in_one <- function(xD){
  D <- f(xD)
  xP <- seq(0.001,0.999,len=201)
  #adjusted variance
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))
  adjusted_var <- em_out[,"VarD_f(x)"]
  plot(xP,em_out[,"VarD_f(x)"],ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=paste("Emulator Adjusted Variance: x_D =",list(xD)))
  #emulator graph
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))  
  plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle=paste("Emulator Output: x_D =",list(xD)))
  #maximum Adjusted Variance
  max_VarD <- max(adjusted_var)
  #mean Adjusted Variance
  mean_VarD <- mean(adjusted_var)
  return(c("max Var_D"=max_VarD,"mean Var_D"=mean_VarD))
}

f <- function(x) cos(3*pi*x)
xD <- c(0,0.1,0.2,0.3,0.6,0.8)
#png(file="poordesign_%d.png",width=2400,height=1600,res=200)
all_in_one(xD)
#dev.off()

#png(file="BIGGER_poordesign_%d.png",width=2100,height=1400,res=200)
all_in_one(xD)
#dev.off()






# Figure 12 ---------------------------------------------------------------

f <- function(x) cos(3*pi*x)

#Investigating how design of runs can affect max and mean emulator variance
#Note that looking at the variance graphs for each:
#Best maximum variance: All variance peaks are roughly equal
#Best mean variance: Variance at the middle is super low and endpoints go up.

finding_best_max_xD <- function(xD){
  D <- f(xD)
  xP <- seq(0.001,0.999,len=201)
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))
  maxvar <- max(em_out[,"VarD_f(x)"])
  return(maxvar)
}

finding_best_mean_xD <- function(xD){
  D <- f(xD)
  xP <- seq(0.001,0.999,len=201)
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))
  meanvar <- mean(em_out[,"VarD_f(x)"])
  return(meanvar)
}

best_max_design <- optim(par=seq(from=0,to=1,length.out=6), fn=finding_best_max_xD, lower=0, upper=1) #equally spaced starting values
max_xD <- best_max_design$par

best_mean_design <- optim(par=seq(from=0,to=1,length.out=6), fn=finding_best_mean_xD, lower=0, upper=1) #equally spaced starting values
mean_xD <- best_mean_design$par

max_to_2dp_all_in_one <- function(xD){
  D <- f(xD)
  xD_to_2dp <- round(xD, digits = 2)
  xP <- seq(0.001,0.999,len=201)
  #adjusted variance
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))
  adjusted_var <- em_out[,"VarD_f(x)"]
  plot(xP,em_out[,"VarD_f(x)"],ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=paste("Minimising Maximum Variance\nEmulator Adjusted Variance: x_D =",list(xD_to_2dp)))
  #emulator graph
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))  
  plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle=paste("Minimising Maximum Variance\nEmulator Output: x_D =",list(xD_to_2dp)))
  #maximum Adjusted Variance
  max_VarD <- max(adjusted_var)
  #mean Adjusted Variance
  mean_VarD <- mean(adjusted_var)
  return(c("max Var_D"=max_VarD,"mean Var_D"=mean_VarD))
}

#png(file="BIGGER_bestmaxvar1D_%d.png",width=2100,height=1400,res=200)
max_to_2dp_all_in_one(max_xD)
#dev.off()

mean_to_2dp_all_in_one <- function(xD){
  D <- f(xD)
  xD_to_2dp <- round(xD, digits = 2)
  xP <- seq(0.001,0.999,len=201)
  #adjusted variance
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))
  adjusted_var <- em_out[,"VarD_f(x)"]
  plot(xP,em_out[,"VarD_f(x)"],ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=paste("Minimising Mean Variance\nEmulator Adjusted Variance: x_D =",list(xD_to_2dp)))
  #emulator graph
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))  
  plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle=paste("Minimising Mean Variance\nEmulator Output: x_D =",list(xD_to_2dp)))
  #maximum Adjusted Variance
  max_VarD <- max(adjusted_var)
  #mean Adjusted Variance
  mean_VarD <- mean(adjusted_var)
  return(c("max Var_D"=max_VarD,"mean Var_D"=mean_VarD))
}

#png(file="BIGGER_bestmeanvar1D_%d.png",width=2100,height=1400,res=200)
mean_to_2dp_all_in_one(mean_xD)
#dev.off()






# Figure 13 ---------------------------------------------------------------

#png(file="var_mean_max_grid.png",width=2400,height=1600,res=200)

#first plot best mean variance
xD <- mean_xD
D <- f(xD)
xP <- seq(0.001,0.999,len=201)
em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))
plot(xP,em_out[,"VarD_f(x)"],ty="l",col="blue",lwd=2.5,
     xlab="Input parameter x",ylab="Output f(x)",main="Adjusted Variance of Emulators with Different Designs")

#add best maximum variance over top
xD <- max_xD
D <- f(xD)
em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))
lines(xP,em_out[,"VarD_f(x)"],col="red",lwd=2.5)

#add variance for equally spaced points over top
xD <- seq(0,1,0.2)
D <- f(xD)
em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))
lines(xP,em_out[,"VarD_f(x)"],col="green",lwd=2.5)
legend('topright',legend=c("Best Mean Variance","Best Maximum Variance","Equally Spaced Runs"),
       lty=c(1,1,1),pch=c(NA,NA,NA),col=c("blue","red","green"),lwd=2.5)


#dev.off()









# Figure 14 ---------------------------------------------------------------

to_2dp_third_new_all_in_one <- function(xD,theta_val,sigma_val){
  D <- f(xD)
  xD_to_2dp <- round(xD, digits = 2)
  xP <- seq(0.001,0.999,len=201)
  #adjusted variance
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=0))
  adjusted_var <- em_out[,"VarD_f(x)"]
  plot(xP,em_out[,"VarD_f(x)"],ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=paste("Theta = ",theta_val,"\nEmulator Variance: x_D =",list(xD_to_2dp)))
  #emulator graph
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=0))  
  plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle=paste("Theta = ",theta_val,"\nEmulator Output: x_D =",list(xD_to_2dp)))
  #maximum Adjusted Variance
  max_VarD <- max(adjusted_var)
  #mean Adjusted Variance
  mean_VarD <- mean(adjusted_var)
  return(c("max Var_D"=max_VarD,"mean Var_D"=mean_VarD))
}
finding_best_max_xD_theta <- function(xD,theta_val){
  D <- f(xD)
  xP <- seq(0.001,0.999,len=201)
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=theta_val,sigma=0.5,E_f=0))
  maxvar <- max(em_out[,"VarD_f(x)"])
  return(maxvar)
}

#theta=0.25
theta.25 <- optim(par=seq(from=0,to=1,length.out=6), fn=finding_best_max_xD_theta, theta_val=0.25, lower=0, upper=1) #equally spaced starting values
theta.25_xD <- theta.25$par

#theta=0.1
theta.1 <- optim(par=seq(from=0,to=1,length.out=6), fn=finding_best_max_xD_theta, theta_val=0.1, lower=0, upper=1) #equally spaced starting values
theta.1_xD <- theta.1$par

#theta=0.35
theta.35 <- optim(par=seq(from=0,to=1,length.out=6), fn=finding_best_max_xD_theta, theta_val=0.35, lower=0, upper=1) #equally spaced starting values
theta.35_xD <- c(0.02737662, 0.18598843, 0.39271606, 0.60728656, 0.81402627, 0.97262147)


#png(file="theta_sensitive.png",width=2400,height=2400,res=200)

par(mfrow=c(3,2))
to_2dp_third_new_all_in_one(theta.1_xD,0.1,0.5)
to_2dp_third_new_all_in_one(theta.25_xD,0.25,0.5)
to_2dp_third_new_all_in_one(theta.35_xD,0.35,0.5)
par(mfrow=c(1,1))

#dev.off()







# Run this before Figure 15 and 16 ----------------------------------------

#load contour plot colour schemes
library(viridisLite)


#needed for efficient emulator function
library(pdist)


f <- function(x) 1.2*cos(2*pi*(x[,1]+0.3)) - 0.7*cos(2*pi*(x[,1]+0.3)*(x[,2]-0.2))

#efficient with nugget: 2D input emulator function 
efficient_nugget_simple_BL_emulator <- function(xP,             # the set of emulator prediction points
                                                xD,             # the run input locations xD
                                                D,              # the run outputs D = (f(x^1),...,f(x^n))
                                                theta = 1,      # the correlation lengths (can be a vector)
                                                sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                                E_f = 0,         # prior expectation of f: E(f(x)) = 0 
                                                using_pdist = 1  # if you have installed pdist package
){
  
  # store length of runs D and number of prediction points xP
  n <- length(D)
  nP <- nrow(xP)       # XXX New V3
  
  # # Rescale each input by theta. Works for different theta for each input and for same theta
  xP <- t(t(xP)/theta)     # XXX New V3: solution to Exercise 8.3, CP 3.
  xD <- t(t(xD)/theta)     # XXX New V3: solution to Exercise 8.3, CP 3.
  
  ### Define Cov structure of f(x): Cov[f(x),f(xdash)], now to act on matrix of distances ###
  # Cov_fx_fxdash <- function(x,xdash) sig^2 * exp(-sum((x-xdash)^2)/theta^2) # XXX Old V2
  Cov_fx_fxdash <- function(dist_matrix) sigma^2 * exp(-(dist_matrix)^2) # XXX New dist V3
  
  
  ### Define 5 objects needed for BL adjustment ###
  # Create E[D] vector
  E_D <- rep(E_f,n)
  
  # Create Var_D matrix:
  # Var_D <- matrix(0,nrow=n,ncol=n)                                        # XXX Old V2
  # for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i,],xD[j,])  # XXX Old V2
  Var_D <- Cov_fx_fxdash( as.matrix(dist(xD)) )                       # XXX New dist V3
  
  
  
  #!!!NEW PART!!!
  #diagonal <- diag(10^-6,n,n) #1 millionth of what we are doing but well above numerical error (nice safe amount)
  diagonal <- diag(0.1^2,n,n) #0.1 as well each green pt quite a bit more uncertain
  Var_D <- Var_D + sigma^2*diagonal
  #END OF NEW PART
  #0.1 or 0.2 closer to 10% of mean
  
  
  
  # Create E[f(x)]
  E_fx <- rep(E_f,nP)
  
  # Create Var_f(x) 
  Var_fx <- rep(sigma^2,nP)
  
  # Create Cov_fx_D row vector now using pdist() function if available, if not use dist()
  # Cov_fx_D <- matrix(0,nrow=1,ncol=n)                       # XXX Old V2
  # for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j,])    # XXX Old V2
  if(using_pdist)  Cov_fx_D <- Cov_fx_fxdash( as.matrix(pdist(xP,xD)) )   # XXX New V3
  if(!using_pdist) 
    Cov_fx_D <- Cov_fx_fxdash( as.matrix(dist(rbind(xP,xD)))[1:nP,(nP+1):(nP+n)] )# XXX NewV3
  
  # find inverse of Var_D using Cholesky decomposition (check Wikipedia if interested!)
  Var_D_inv <- chol2inv(chol(Var_D))     # more efficient way to calculate inverse of Cov mat
  
  ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
  cov_fx_D_Var_D_inv <- Cov_fx_D %*% Var_D_inv  # Need this twice so pre-calculate here
  ED_fx   <-  E_fx + cov_fx_D_Var_D_inv %*% (D - E_D) # adj. expectation of ALL xP pts at once
  # VarD_fx <-  Var_fx - cov_fx_D_Var_D_inv %*% t(Cov_fx_D)       # XXX Old V2
  VarD_fx   <-  Var_fx - apply(cov_fx_D_Var_D_inv * Cov_fx_D,1,sum) # fast way to get diagonals 
  # and hence all nP variances (Note: does not do full np x np covariance matrix)
  
  ### return emulator expectation and variance ###
  return(cbind("ExpD_f(x)"=c(ED_fx),"VarD_f(x)"=VarD_fx))  
  
}


#produces adjusted variance graph
efficient_nugget_vector_input_var_graph_only <- function(xD_vec, 
                                                         grid_length=50, 
                                                         theta_val=0.5, 
                                                         sigma_val=0.5, 
                                                         E_f_val=0)
{
  par(mar=c(4,4,3,3))
  var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
  
  xD <- matrix(xD_vec,ncol=2,byrow=TRUE)
  D <- f(xD)
  
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  em_out <- efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val)   
  
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")
  
}








#creating contour plot
emul_fill_cont <- function(cont_mat,            # matrix of values we want contour plot of 
                           cont_levs=NULL,      # contour levels (NULL: automatic selection)
                           nlev=20,             # approx no. of contour levels for auto select  
                           plot_xD=TRUE,        # plot the design runs TRUE or FALSE
                           xD=NULL,             # the design points if needed
                           xD_col="green",      # colour of design runs
                           x_grid,              # grid edge locations that define xP
                           ...                  # extra arguments passed to filled.contour
){
  
  ### Define contour levels if necessary ###
  if(is.null(cont_levs)) cont_levs <- pretty(cont_mat,n=nlev)     
  
  ### create the filled contour plot ###
  filled.contour(x_grid,x_grid,cont_mat,levels=cont_levs,xlab="x1",ylab="x2",...,  
                 plot.axes={axis(1);axis(2)                 # sets up plotting in contour box
                   contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs,lwd=0.8)   # plot contour lines
                   if(plot_xD) points(xD,pch=21,col=1,bg=xD_col,cex=1.5)})  # plot design points 
  #if we want design pts ie plot_xD=TRUE
}


#plotting contour plots for 
#     adjusted expectation
#     adjusted variance 

efficient_nugget_all_graphs <- function(xD_matrix, 
                                        D_vector,
                                        grid_length=50, 
                                        theta_val=0.5, 
                                        sigma_val=0.5, 
                                        E_f_val=0){
  
  par(mar=c(4,4,3,3))
  
  #colours
  exp_cols <- turbo #see ?viridis
  var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
  
  
  #run locations and their outputs
  xD <- xD_matrix 
  D <- D_vector 
  
  
  #grid needed for emulation
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  
  #emulation
  em_out <- efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val)  
  
  
  #extracting the emulator expectation and variance
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  
  #graphing
  emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")
}





#plotting contour plots using the "efficient_nugget_all_graphs" function
#     this function has a vector input
#     this is useful after using optim

efficient_nugget_vector_input_all_graphs <- function(xD_vec,
                                                     D_vector,
                                                     grid_length=50, 
                                                     theta_val=0.5, 
                                                     sigma_val=0.5, 
                                                     E_f_val=0)
{
  par(mar=c(4,4,3,3))
  xD <- matrix(xD_vec,ncol=2,byrow=TRUE)
  efficient_nugget_all_graphs(xD_matrix = xD, D_vector, grid_length, theta_val, sigma_val, E_f_val)
}



#optimiser function for min quantile
efficient_nugget_2D_finding_best_quantile_xD <- function(xD_vector,
                                                         grid_length=50,
                                                         theta_val=0.5,
                                                         sigma_val=0.5,
                                                         E_f_val=0,
                                                         quantile_we_choose=0.9){ #will optimise first argument
  xD <- matrix(xD_vector,ncol=2,byrow=TRUE)
  D <- f(xD)
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  em_out <- efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val)
  quantile <- quantile(em_out[,"VarD_f(x)"], probs=quantile_we_choose)
  return(quantile)
}






#optimiser function for min max
efficient_nugget_2D_finding_best_max_xD <- function(xD_vector,
                                                    grid_length=50,
                                                    theta_val=0.5,
                                                    sigma_val=0.5,
                                                    E_f_val=0){ #will optimise first argument
  xD <- matrix(xD_vector,ncol=2,byrow=TRUE)
  D <- f(xD)
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  em_out <- efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val)
  maxvar <- max(em_out[,"VarD_f(x)"])
  return(maxvar)
}





#optimiser function for min mean
efficient_nugget_2D_finding_best_mean_xD <- function(xD_vector,
                                                     grid_length=50,
                                                     theta_val=0.5,
                                                     sigma_val=0.5,
                                                     E_f_val=0){ #will optimise first argument
  xD <- matrix(xD_vector,ncol=2,byrow=TRUE)
  D <- f(xD)
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  em_out <- efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val)
  meanvar <- mean(em_out[,"VarD_f(x)"])
  return(meanvar)
}






# Figure 15 ---------------------------------------------------------------

starting_pts_6 <- c(0.1,0.2,
                    0.5,0.3,
                    0.9,0.1,
                    0.15,0.9,
                    0.4,0.7,
                    0.8,0.9)
#just to check nicely spread
par(mfrow=c(1,1))
plot(matrix(starting_pts_6,ncol=2,byrow = TRUE),xlim=c(0,1),ylim=c(0,1))

#mean
mean_6 <- optim(method = "L-BFGS-B",
                
                par = starting_pts_6,
                
                fn = efficient_nugget_2D_finding_best_mean_xD,
                
                
                grid_length=10,
                theta_val=0.5,
                sigma_val=1,
                E_f_val=0,
                
                
                lower=0,
                upper=1)

pts_mean_6 <- mean_6$par

#png(file="week6_6_points_mean.png",width=1400,height=1200,res=200)
efficient_nugget_vector_input_var_graph_only(pts_mean_6)
#dev.off()







# Figure 16 ---------------------------------------------------------------

par(mfrow=c(1,1))

#roughly equally spaced starting points (not exact grid so optim does not fall into any symmetric traps)
starting_pts_9 <- c(0.1,0.15,
                    0.5,0.2,
                    0.8,0.05,
                    0.2,0.5,
                    0.55,0.6,
                    0.9,0.45,
                    0.15,0.9,
                    0.45,0.85,
                    0.8,0.95)
#just to check nicely spread
par(mfrow=c(1,1))
#png(file="week6_9_points_starting.png",width=1100,height=1200,res=200)
plot(matrix(starting_pts_9,ncol=2,byrow = TRUE),xlim=c(0,1),ylim=c(0,1),
     xlab="x1",ylab="x2",main="Initial Values")
#dev.off()


#max
max_9 <- optim(method = "L-BFGS-B",
               
               par = starting_pts_9,
               
               fn = efficient_nugget_2D_finding_best_max_xD,
               
               
               grid_length=10,
               theta_val=0.5,
               sigma_val=1,
               E_f_val=0,
               
               
               lower=0,
               upper=1)

pts_max_9 <- max_9$par

#png(file="week6_9_points_max.png",width=1400,height=1200,res=200)
efficient_nugget_vector_input_var_graph_only(pts_max_9)
#dev.off()

#mean
mean_9 <- optim(method = "L-BFGS-B",
                
                par = starting_pts_9,
                
                fn = efficient_nugget_2D_finding_best_mean_xD,
                
                
                grid_length=10,
                theta_val=0.5,
                sigma_val=1,
                E_f_val=0,
                
                
                lower=0,
                upper=1)

pts_mean_9 <- mean_9$par

#png(file="week6_9_points_mean.png",width=1400,height=1200,res=200)
efficient_nugget_vector_input_var_graph_only(pts_mean_9)
#dev.off()

#90 quantile
quantile90_9 <- optim(method = "L-BFGS-B",
                      
                      par = starting_pts_9,
                      
                      fn = efficient_nugget_2D_finding_best_quantile_xD,
                      
                      
                      grid_length=60,
                      theta_val=0.5,
                      sigma_val=1,
                      E_f_val=0,
                      
                      quantile_we_choose=0.7,
                      
                      
                      lower=0,
                      upper=1)

pts_quantile90_9 <- quantile90_9$par

#png(file="week6_9_points_quantile.png",width=1400,height=1200,res=200)
efficient_nugget_vector_input_var_graph_only(pts_quantile90_9)
#dev.off()







# End of Section 3 --------------------------------------------------------
