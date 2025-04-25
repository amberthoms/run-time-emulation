
# Appendices --------------------------------------------------------------

# Figures 35-43

# Functions ---------------------------------------------------------------

# Figure 35 ---------------------------------------------------------------
## nailing down a function -------------------------------------------------
NEW_VAR_D_simple_BL_emulator_v1 <- function(x,              # the emulator prediction point
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
  
  #!!!NEW PART!!!
  diagonal <- diag(10^-6,n,n)
  Var_D <- Var_D + sigma^2*diagonal
  #END OF NEW PART
  
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



naildown_plot_BL_emulator_V1 <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))^T
    maintitle=NULL  # title of the plot, blank if NULL
){
  
  ### plot emulator output ###
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.4,1.4),ty="l",col="blue",lwd=2.5, lty=1,
       xlab="Input parameter x",ylab="Output f(x)",main=maintitle) # main argument: plot title (INPUT TO SPECIFY)
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  
  ### plot true function: we would not normally be able to do this! ###
  lines(xP,f(xP),lwd=2,lty=2)
  
  ### Plot the runs ### 
  points(xD,D,pch=21,col=1,bg="green",cex=1.5) #PLOT THE TRUE RUNS OF THE EXPENSIVE MODEL
  legend('topright',legend=c("Emulator Expectation","Emulator Prediction Interval",
                             "True function f(x)","Model Evaluations"),
         lty=c(1,1,2,NA),pch=c(NA,NA,NA,16),col=c("blue","red",1,"green"),lwd=2.5,pt.cex=1.3)
}


xD <- seq(0,1,len=30) 
#LUXURY TO HAVE THIS MANY PTS
#MAYBE HUNDREDS OF THOUSANDS NEEDED TO NAIL HIGHER DIM FN

D <- f(xD)
xP <- seq(0.001,0.999,len=201)

em_out <- t(sapply(xP,NEW_VAR_D_simple_BL_emulator_v1,xD=xD,D=D,theta=0.05,sigma=0.5,E_f=0))

#png(file="naildownfunction.png",width=2400,height=1900,res=200)

naildown_plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Emulator Output: Nailing Down the Function")

#dev.off()






# Figure 36 ---------------------------------------------------------------

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


##############################################################################################.
### Function to plot simple emulator output
plot_BL_emulator_V1 <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))^T
    maintitle=NULL  # title of the plot, blank if NULL
){
  
  ### plot emulator output ###
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.2,1.2),ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=maintitle) # main argument: plot title (INPUT TO SPECIFY)
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="red",lwd=2.5)
  
  ### plot true function: we would not normally be able to do this! ###
  lines(xP,f(xP),lwd=2,lty=1)
  
  ### Plot the runs ### 
  points(xD,D,pch=21,col=1,bg="green",cex=1.5) #PLOT THE TRUE RUNS OF THE EXPENSIVE MODEL
  legend('topright',legend=c("Emulator Expectation","Emulator Prediction Interval",
                             "True function f(x)","Model Evaluations"),
         lty=c(1,1,1,NA),pch=c(NA,NA,NA,16),col=c("blue","red",1,"green"),lwd=2.5,pt.cex=1.3)
}
### Function to plot simple emulator output
##############################################################################################.


### Define actual computer model/simulator ###
f <- function(x) cos(3*pi*x) #function named f with one input x
#f <- function(x) sin(2*pi*x)



## design sensitive to sigma? ----------------------------------------------
#No.

finding_best_max_xD_theta_sigma <- function(xD,theta_val,sigma_val){
  D <- f(xD)
  xP <- seq(0.001,0.999,len=201)
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=0))
  maxvar <- max(em_out[,"VarD_f(x)"])
  return(maxvar)
}

#sigma 0.5,1,0.2
sigma.5 <- optim(par=seq(from=0,to=1,length.out=6), fn=finding_best_max_xD_theta_sigma, theta_val=0.25, sigma_val=0.5, lower=0, upper=1) #equally spaced starting values
sigma1 <- optim(par=seq(from=0,to=1,length.out=6), fn=finding_best_max_xD_theta_sigma, theta_val=0.25, sigma_val=1, lower=0, upper=1)
sigma.2 <- optim(par=seq(from=0,to=1,length.out=6), fn=finding_best_max_xD_theta_sigma, theta_val=0.25, sigma_val=0.2, lower=0, upper=1)

sigma.2$par
sigma.5$par
sigma1$par

#these are very, very similar - 
#     design independent of sigma - just scales up and down

to_2dp_second_new_all_in_one <- function(xD,theta_val,sigma_val){
  D <- f(xD)
  xD_to_2dp <- round(xD, digits = 2)
  xP <- seq(0.001,0.999,len=201)
  #adjusted variance
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=0))
  adjusted_var <- em_out[,"VarD_f(x)"]
  plot(xP,em_out[,"VarD_f(x)"],ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=paste("Sigma = ",sigma_val,"\nEmulator Variance: x_D =",list(xD_to_2dp)))
  #emulator graph
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=0))  
  plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle=paste("Sigma = ",sigma_val,"\nEmulator Output: x_D =",list(xD_to_2dp)))
  #maximum Adjusted Variance
  max_VarD <- max(adjusted_var)
  #mean Adjusted Variance
  mean_VarD <- mean(adjusted_var)
  return(c("max Var_D"=max_VarD,"mean Var_D"=mean_VarD))
}

#the following graphs show that sigma just scales emulator variance

#png(file="sigma_sensitive.png",width=2400,height=2400,res=200)

par(mfrow=c(3,2))
to_2dp_second_new_all_in_one(sigma.2$par,0.25,0.2)
to_2dp_second_new_all_in_one(sigma.5$par,0.25,0.5)
to_2dp_second_new_all_in_one(sigma1$par,0.25,1)
par(mfrow=c(1,1))

#dev.off()






# Figure 37 ---------------------------------------------------------------


## run before alpha --------------------------------------------------------




rm(list = ls())

library(viridisLite)
library(pdist)

nugget_argument_efficient_nugget_all_graphs <- function(xD_matrix, 
                                                        D_vector,
                                                        grid_length=50, 
                                                        theta_val=0.5, 
                                                        sigma_val=0.5, 
                                                        E_f_val=0,
                                                        
                                                        #new part:
                                                        nugget_val=10^(-6)){
  
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
  #new part (adding nugget as argument):
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)  
  
  
  #extracting the emulator expectation and variance
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  
  #graphing
  emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")
}





nugget_argument_efficient_nugget_simple_BL_emulator <- function(xP,             # the set of emulator prediction points
                                                                xD,             # the run input locations xD
                                                                D,              # the run outputs D = (f(x^1),...,f(x^n))
                                                                theta = 1,      # the correlation lengths (can be a vector)
                                                                sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                                                E_f = 0,        # prior expectation of f: E(f(x)) = 0 
                                                                using_pdist = 1,# if you have installed pdist package
                                                                
                                                                #new part:
                                                                nugget = 10^(-6) #must input nugget 0.1^2 
                                                                #squared because working with sigma squared
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
  
  
  
  
  
  
  #deterministic: manage numerical error
  #diagonal <- diag(10^-6,n,n) #1 millionth of what we are doing but well above numerical error (nice safe amount)
  #stochastic: manage stochasticity
  #diagonal <- diag(0.1^2,n,n) #each green pt quite a bit more uncertain
  
  #new part:
  diagonal <- diag(nugget,n,n)
  Var_D <- Var_D + sigma^2*diagonal
  
  #note that nugget 0.1^2 or 0.2^2 closer to 10% of mean of D (which is roughly what we want)
  #could put in directly but will have to change by hand every time change emulator parameters
  
  
  
  
  
  
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


dataframe_sim_run_times <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/simulation_overview.txt")
matrix_sim_run_times <- as.matrix(dataframe_sim_run_times)
colnames(matrix_sim_run_times) <- c("Emulator Coordinate x", 
                                    "Emulator Coordinate y", 
                                    "Runtime [CPU-Hours]", 
                                    "Batch Number", 
                                    "fmin",	
                                    "fmax",	
                                    "rho0 [cm^-3]",
                                    "T_AGN [K]")
x <- matrix_sim_run_times[,"Emulator Coordinate x"]
y <- matrix_sim_run_times[,"Emulator Coordinate y"]

xD <- matrix(data=c(x,y),ncol=2,byrow = FALSE)
D <- log10(matrix_sim_run_times[,"Runtime [CPU-Hours]"])

emulate_log_run_time_Tc <- function(new_batch_vector,
                                    xD_matrix,
                                    D_vector,
                                    
                                    theta_val,
                                    sigma_val,
                                    E_f_val,
                                    
                                    nugget_val=10^(-6))
  
{
  xD <- xD_matrix
  D <- D_vector
  new_batch_matrix <- matrix(data=new_batch_vector,ncol=2,byrow = TRUE)
  xP <- new_batch_matrix
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)
  E_D_fx_mat <- em_out[,"ExpD_f(x)"] 
  normal_times <- 10^(E_D_fx_mat)
  return(sum(normal_times))
}

var_graph <- function(xD_vec_original,
                      xD_vec_new,
                      grid_length=50, 
                      theta_val=0.5, 
                      sigma_val=0.5, 
                      E_f_val=0,
                      
                      nugget_val=10^(-6))
{
  par(mar=c(4,4,3,3))
  var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
  
  xD_vec <- c(xD_vec_original,xD_vec_new)
  xD <- matrix(xD_vec,ncol=2,byrow=TRUE)
  D <- f(xD)
  
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)   
  
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]",
                 #new part:
                 xD_col=rep(c("purple","pink"),c(length(xD_vec_original)/2,length(xD_vec_new)/2)) #/2 because x and y coords in vector so number of points is length/2
  )
  
}

#create random function (won't affect variance but function needs it)
f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))


library(MASS)

noise_function <- function(batch_size,
                           number_ICs,
                           sigma_noise){
  var_matrix <- matrix(c(sigma_noise^2,0,
                         0,sigma_noise^2),
                       nrow=2,
                       ncol=2,
                       byrow=TRUE
  )
  noise2D_matrix <- mvrnorm(n=batch_size*(number_ICs-1),
                            mu=rep(0,2),
                            Sigma=var_matrix)
  noise2D_vector <- c(rep(0,batch_size*2),c(t(noise2D_matrix)))
  return(noise2D_vector)
}

new_optim_fixed_runs <- function(newpts, #will optimise first argument
                                 oldpts,
                                 grid_length=50,
                                 theta_val2=0.5,
                                 sigma_val2=0.5,
                                 E_f_val2=0,
                                 nugget_val2=10^(-6),
                                 true_runs_vec,
                                 alpha=1,
                                 
                                 runtime_emulator_theta,
                                 runtime_emulator_sigma,
                                 runtime_emulator_E,
                                 runtime_emulator_nugget
){ 
  oldpts_matrix <- matrix(oldpts,ncol=2,byrow=TRUE) 
  newpts_matrix <- matrix(newpts,ncol=2,byrow=TRUE)
  xD <- rbind(oldpts_matrix,newpts_matrix)
  D <- f(xD)
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val2,sigma=sigma_val2,E_f=E_f_val2,nugget=nugget_val2)
  v90 <- quantile(em_out[,"VarD_f(x)"],0.9)
  Tc <- emulate_log_run_time_Tc(new_batch_vector=newpts,
                                xD_matrix=oldpts_matrix,
                                D_vector=true_runs_vec,
                                theta_val=runtime_emulator_theta,
                                sigma_val=runtime_emulator_sigma,
                                E_f_val=runtime_emulator_E,
                                nugget_val=runtime_emulator_nugget)
  final <- (1/2)*log10(v90) + alpha*log10(Tc)
  return(final)
  
}

new_var_grid_graph <- function(want_var_graph,
                               want_grid,
                               want_graph,
                               want_best_var_graphs,
                               
                               batch0,
                               batch0_evaluated, #evaluated run times
                               
                               grid_length,
                               
                               batch_size, #4, 6, 8 or 10
                               
                               theta_emulator,
                               sigma_emulator,
                               E_emulator,
                               nugget_emulator, #input as squared
                               
                               run_time_emulator_theta,
                               run_time_emulator_sigma,
                               run_time_emulator_E,
                               run_time_emulator_nugget,
                               
                               alphas, #vector of alphas we want
                               percentile, #percentile of variance used in optim 
                               number_ICs, #how many different starting points we want
                               sigma_noise #N(0,sigma_noise^2) for random noise
){
  
  
  
  #determine initial starting points 
  #dependent on batch_size (number of new points) we want
  if (batch_size==4){
    starting_pts <- c(0.9,0.85,
                      0.6,0.95, 
                      0.4,0.05,
                      0.1,0.35)
  }
  if (batch_size==6){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05)
  }
  if (batch_size==8){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6)
  }
  if (batch_size==10){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6,
                      
                      0.2,0.8,
                      0.8,0.1)
  }
  
  
  
  
  #this creates "number_ICs" lots of initial starting points
  #1st is with no noise
  #2nd onwards are with 2D multivariate normal noise mean (0,0) var=(sigma_noise^2 , 0
  #                                                                  0 , sigma_noise^2)
  noise_vector <- noise_function(batch_size,
                                 number_ICs,
                                 sigma_noise)
  
  starting_pts_without_ICs <- rep(starting_pts,number_ICs)
  starting_pts_with_ICs <- starting_pts_without_ICs+noise_vector
  
  
  
  
  #this ensures any noisy initial values below 0 are changed to 0.01
  #similarly above 1 changed to 0.99
  starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs<0,0.01)
  starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs>1,0.99)
  
  
  
  
  #now we find the optim values for each alpha for each IC
  #we save the parameters and the optimised value
  optim_par <- rep(0,length(alphas)*number_ICs*batch_size*2)
  optim_value <- rep(0,length(alphas)*number_ICs)
  z <- 1
  for (i in 1:length(alphas)){
    for (j in 1:number_ICs){
      optim <- optim(method = "L-BFGS-B",
                     
                     par=starting_pts_with_ICs[(2*(j-1)*batch_size+1):(2*j*batch_size)], 
                     
                     fn=new_optim_fixed_runs,
                     
                     oldpts=batch0,
                     grid_length=grid_length,
                     theta_val2=theta_emulator,
                     sigma_val2=sigma_emulator,
                     E_f_val2=E_emulator,
                     nugget_val2=nugget_emulator,
                     true_runs_vec=batch0_evaluated,
                     alpha=alphas[i],
                     
                     runtime_emulator_theta=run_time_emulator_theta,
                     runtime_emulator_sigma=run_time_emulator_sigma,
                     runtime_emulator_E=run_time_emulator_E,
                     runtime_emulator_nugget=run_time_emulator_nugget,
                     
                     lower=0,
                     upper=1)
      optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)] <- optim$par
      optim_value[z] <- optim$value
      z <- z+1
    }
  }
  
  
  
  
  #this produces var graph showing locations of optimised points for every IC for every alpha
  #warning - this can produce hundreds of graphs!
  if (want_var_graph==TRUE){
    for (z in 1:(length(alphas)*number_ICs)){
      var_graph(xD_vec_original=batch0,
                xD_vec_new=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                grid_length=grid_length, 
                theta_val=theta_emulator, 
                sigma_val=sigma_emulator, 
                E_f_val=E_emulator,
                nugget_val=nugget_emulator)
    }
  }
  
  
  
  
  #we calculate v90 and Tc for each set of optimised locations
  if (want_grid==TRUE | want_graph==TRUE){
    #calculate V90 and Tc for each set of locations
    oldpts_matrix <- matrix(batch0,ncol=2,byrow=TRUE) 
    v90 <- rep(0,length(alphas)*number_ICs)
    Tc <- rep(0,length(alphas)*number_ICs)
    for (z in 1:(length(alphas)*number_ICs)){ #look at this section again... particularly the v90
      newpts_matrix <- matrix(optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],ncol=2,byrow=TRUE)
      xD <- rbind(oldpts_matrix,newpts_matrix)
      D <- f(xD)
      x_grid <- seq(-0.001,1.001,len=grid_length)
      xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
      em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_emulator,sigma=sigma_emulator,E_f=E_emulator,nugget=nugget_emulator)
      v90[z] <- quantile(em_out[,"VarD_f(x)"],percentile)
      Tc[z] <- emulate_log_run_time_Tc(new_batch_vector=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                                       xD_matrix=oldpts_matrix,
                                       D_vector=batch0_evaluated,
                                       theta_val=run_time_emulator_theta,
                                       sigma_val=run_time_emulator_sigma,
                                       E_f_val=run_time_emulator_E,
                                       nugget_val=run_time_emulator_nugget) #take care this should use the log emulator characterisitics!!! May have to change in future
    }
    
    
  }
  
  
  
  
  #for each alpha we find best locations (over the ICs) based on which gives the lowest optim value
  minimums_indices <- rep(0,length(alphas))
  for (k in 1:length(alphas)){
    minimums_indices[k] <- (k-1)*number_ICs+which.min(optim_value[((k-1)*number_ICs+1):(k*number_ICs)])
  }
  #pars_for_grid <- list_optim_par[[minimums_indices]] #this of length length(alphas)
  v90_for_grid <- v90[minimums_indices]
  Tc_for_grid <- Tc[minimums_indices]
  
  
  
  
  #a grid 
  #for each row e.g. 1/2 evaluate the relative criterion log10U = 1/2*log10V_{90}+1/2*log10T_{c} for the optim design using each alpha
  #if optim is performing/behaving correctly the diagonals should be the minimum of each row 
  if (want_grid==TRUE){
    #minimum for each different starting conditions
    matrix <- matrix(0,ncol=length(alphas),nrow=length(alphas))
    for (i in 1:length(alphas)){
      for (j in 1:length(alphas)){
        matrix[i,j] <- 1/2*log10(v90_for_grid[j])+alphas[i]*log10(Tc_for_grid[j])
      }
    }
    colnames(matrix) <- alphas
    rownames(matrix) <- alphas
    print(as.table(matrix))
  }
  
  
  
  
  #plot Tc on x axis against V90 on y axis for each alpha
  if (want_graph==TRUE){
    #graph
    par(mfrow=c(1,1))
    plot(Tc,v90,
         xlab="Tc",
         ylab="V90")
    #label the plot with fewer points
    plot(Tc_for_grid,v90_for_grid,
         xlab="Tc",
         ylab="V90")
    text(Tc_for_grid,v90_for_grid,label=round(alphas,2),cex= 0.7,pos=3)
  }
  
  
  
  
  #plots the variance graphs (showing where optimised locations are) only for the best locations for each alpha
  #alpha given in bottom left corner
  if (want_best_var_graphs==TRUE){
    for (z in 1:length(alphas)){
      var_graph(xD_vec_original=batch0,
                xD_vec_new=optim_par[((batch_size*2)*(minimums_indices[z]-1)+1):((batch_size*2)*minimums_indices[z])],
                grid_length=grid_length, 
                theta_val=theta_emulator, 
                sigma_val=sigma_emulator, 
                E_f_val=E_emulator,
                nugget_val=nugget_emulator)
      title(paste("Alpha =",round(alphas[z],2)),adj=0,line=-26)
    }
  }
  
  
  
  
  #return some things we may want to look at - especially the table (easy to save this way)
  return(list("Tc"=Tc,"v90"=v90,"Tc_for_grid"=Tc_for_grid,"v90_for_grid"=v90_for_grid,"minimums_indices"=minimums_indices,"optim_par"=optim_par,"optim_value"=optim_value,"table"=as.table(matrix)))
}



## end of run before alpha. ------------------------------------------------



## alpha: Tc against V90, all ICs, best designs ----------------------------
#6 points
#11 alphas
#10 ICs each
#different alphas
#theta=0.5 or 0.3...

#png(file="q_alphas_ICs_%d.png",width=1400,height=1200,res=200)
q_alphas_ICs <- new_var_grid_graph(want_var_graph=FALSE,
                                   want_grid=TRUE,
                                   want_graph=TRUE,
                                   want_best_var_graphs=TRUE,
                                   
                                   batch0=c(t(xD)),
                                   batch0_evaluated=D,
                                   
                                   grid_length=60,
                                   
                                   batch_size=6, #4, 6, 8 or 10
                                   
                                   theta_emulator=1/2,
                                   sigma_emulator=sd(D),
                                   E_emulator=mean(D),
                                   nugget_emulator=0.001^2,
                                   
                                   run_time_emulator_theta=1,
                                   run_time_emulator_sigma=sd(D),
                                   run_time_emulator_E=mean(D),
                                   run_time_emulator_nugget=0.145064^2,
                                   
                                   alphas=c(0,0.1,1/3,1/2,2/3,0.9,1,1.5,2,3,4), #vector of alphas we want
                                   percentile=0.9, #percentile of variance used in optim 
                                   number_ICs=10, #how many different starting points we want
                                   sigma_noise=0.1)

#dev.off()



# Run this before Figure 38 and Figure 39 ---------------------------------

## run before ceiling... ---------------------------------------------------
rm(list = ls())


library(viridisLite)
library(pdist)

nugget_argument_efficient_nugget_all_graphs <- function(xD_matrix, 
                                                        D_vector,
                                                        grid_length=50, 
                                                        theta_val=0.5, 
                                                        sigma_val=0.5, 
                                                        E_f_val=0,
                                                        
                                                        #new part:
                                                        nugget_val=10^(-6)){
  
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
  #new part (adding nugget as argument):
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)  
  
  
  #extracting the emulator expectation and variance
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  
  #graphing
  emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")
}





nugget_argument_efficient_nugget_simple_BL_emulator <- function(xP,             # the set of emulator prediction points
                                                                xD,             # the run input locations xD
                                                                D,              # the run outputs D = (f(x^1),...,f(x^n))
                                                                theta = 1,      # the correlation lengths (can be a vector)
                                                                sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                                                E_f = 0,        # prior expectation of f: E(f(x)) = 0 
                                                                using_pdist = 1,# if you have installed pdist package
                                                                
                                                                #new part:
                                                                nugget = 10^(-6) #must input nugget 0.1^2 
                                                                #squared because working with sigma squared
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
  
  
  
  
  
  
  #deterministic: manage numerical error
  #diagonal <- diag(10^-6,n,n) #1 millionth of what we are doing but well above numerical error (nice safe amount)
  #stochastic: manage stochasticity
  #diagonal <- diag(0.1^2,n,n) #each green pt quite a bit more uncertain
  
  #new part:
  diagonal <- diag(nugget,n,n)
  Var_D <- Var_D + sigma^2*diagonal
  
  #note that nugget 0.1^2 or 0.2^2 closer to 10% of mean of D (which is roughly what we want)
  #could put in directly but will have to change by hand every time change emulator parameters
  
  
  
  
  
  
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



dataframe_sim_run_times <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/simulation_overview.txt")
matrix_sim_run_times <- as.matrix(dataframe_sim_run_times)
colnames(matrix_sim_run_times) <- c("Emulator Coordinate x", 
                                    "Emulator Coordinate y", 
                                    "Runtime [CPU-Hours]", 
                                    "Batch Number", 
                                    "fmin",	
                                    "fmax",	
                                    "rho0 [cm^-3]",
                                    "T_AGN [K]")
x <- matrix_sim_run_times[,"Emulator Coordinate x"]
y <- matrix_sim_run_times[,"Emulator Coordinate y"]
plot(x,y)

xD <- matrix(data=c(x,y),ncol=2,byrow = FALSE)
D <- log10(matrix_sim_run_times[,"Runtime [CPU-Hours]"])

emulate_log_run_time_Tc <- function(new_batch_vector,
                                    xD_matrix,
                                    D_vector,
                                    
                                    theta_val,
                                    sigma_val,
                                    E_f_val,
                                    
                                    nugget_val=10^(-6))
  
{
  xD <- xD_matrix
  D <- D_vector
  new_batch_matrix <- matrix(data=new_batch_vector,ncol=2,byrow = TRUE)
  xP <- new_batch_matrix
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)
  E_D_fx_mat <- em_out[,"ExpD_f(x)"] 
  normal_times <- 10^(E_D_fx_mat)
  return(sum(normal_times))
}


var_graph <- function(xD_vec_original,
                      xD_vec_new,
                      grid_length=50, 
                      theta_val=0.5, 
                      sigma_val=0.5, 
                      E_f_val=0,
                      
                      nugget_val=10^(-6))
{
  par(mar=c(4,4,3,3))
  var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
  
  xD_vec <- c(xD_vec_original,xD_vec_new)
  xD <- matrix(xD_vec,ncol=2,byrow=TRUE)
  D <- f(xD)
  
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)   
  
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]",
                 #new part:
                 xD_col=rep(c("purple","pink"),c(length(xD_vec_original)/2,length(xD_vec_new)/2)) #/2 because x and y coords in vector so number of points is length/2
  )
  
}

#create random function (won't affect variance but function needs it)
f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))





library(MASS)

noise_function <- function(batch_size,
                           number_ICs,
                           sigma_noise){
  var_matrix <- matrix(c(sigma_noise^2,0,
                         0,sigma_noise^2),
                       nrow=2,
                       ncol=2,
                       byrow=TRUE
  )
  noise2D_matrix <- mvrnorm(n=batch_size*(number_ICs-1),
                            mu=rep(0,2),
                            Sigma=var_matrix)
  noise2D_vector <- c(rep(0,batch_size*2),c(t(noise2D_matrix)))
  return(noise2D_vector)
}






sigmoidal_optim_function <- function(newpts, #will optimise first argument
                                     oldpts,
                                     grid_length=50,
                                     theta_val2=0.5,
                                     sigma_val2=0.5,
                                     E_f_val2=0,
                                     nugget_val2=10^(-6),
                                     true_runs_vec,
                                     
                                     ceiling,
                                     freeze,
                                     
                                     runtime_emulator_theta,
                                     runtime_emulator_sigma,
                                     runtime_emulator_E,
                                     runtime_emulator_nugget
){
  oldpts_matrix <- matrix(oldpts,ncol=2,byrow=TRUE) 
  newpts_matrix <- matrix(newpts,ncol=2,byrow=TRUE)
  xD <- rbind(oldpts_matrix,newpts_matrix)
  D <- f(xD)
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val2,sigma=sigma_val2,E_f=E_f_val2,nugget=nugget_val2)
  v90 <- quantile(em_out[,"VarD_f(x)"],0.9)
  Tc <- emulate_log_run_time_Tc(new_batch_vector=newpts,
                                xD_matrix=oldpts_matrix,
                                D_vector=true_runs_vec,
                                theta_val=runtime_emulator_theta,
                                sigma_val=runtime_emulator_sigma,
                                E_f_val=runtime_emulator_E,
                                nugget_val=runtime_emulator_nugget)
  
  if (freeze==FALSE){
    a <- 5 #steepness
    b <- 3 #height
    c <- ceiling #log10(ceiling) x value halfway up the sigmoidal
    d <- 0.5 #steepness of linear term added on 
  }
  
  if (freeze==TRUE){
    a <- 10 #steepness - THIS IS THE ONLY DIFFERENCE
    b <- 3 #height
    c <- ceiling #log10(ceiling) x value halfway up the sigmoidal
    d <- 0.5 #steepness of linear term added on
  }
  
  x <- Tc
  
  if (Tc<=ceiling){
    sigmoidal <- b/(1+exp(-a*(x-c)))
  }
  if (Tc>ceiling){
    sigmoidal <- b/(1+exp(-a*(x-c))) + d*x - d*c
  }
  
  final <- (1/2)*log10(v90) + sigmoidal
  return(final)
}

# #checking sigmoidal part works ok
# sig <- function(log10Tc,
#                 log10ceiling
#                 ){
#   a <- 10 #5 or 10
#   b <- 3
#   c <- log10ceiling
#   d <- 0.5
#   x <- log10Tc
#   if (log10Tc<=log10ceiling){
#     sigmoidal <- b/(1+exp(-a*(x-c)))
#   }
#   if (log10Tc>log10ceiling){
#     sigmoidal <- b/(1+exp(-a*(x-c))) + d*x - d*c
#   }
#   return(sigmoidal)
# }
# x <- seq(0,10,by=0.01)
# which(x==4)
# plot(x,sapply(x,FUN=sig,log10ceiling=7),ylim=c(0,6))
# y <- sapply(x,FUN=sig,log10ceiling=7)
# y[601]

#Function using sigmoidal criteria
everything_function_sigmoidal <- function(want_all_var_graphs,
                                          want_only_best_var_graph,
                                          want_freeze, #initialise design at melt design found now with harsher sigmoidal
                                          
                                          batch0,
                                          batch0_evaluated, #evaluated run times
                                          
                                          grid_length,
                                          
                                          batch_size, #4, 6, 8 or 10 now 5,7 also allowed
                                          
                                          theta_emulator,
                                          sigma_emulator,
                                          E_emulator,
                                          nugget_emulator, #input as squared
                                          
                                          run_time_emulator_theta,
                                          run_time_emulator_sigma,
                                          run_time_emulator_E,
                                          run_time_emulator_nugget,
                                          
                                          percentile, #percentile of variance used in optim 
                                          number_ICs, #how many different starting points we want
                                          sigma_noise, #N(0,sigma_noise^2) for random noise
                                          
                                          ceiling_value #NOT LOGGED max total run time cannot exceed
){
  
  #create random function (won't affect variance but function needs it)
  f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))
  
  
  
  #determine initial starting points 
  #dependent on batch_size (number of new points) we want
  if (batch_size==4){
    starting_pts <- c(0.9,0.85,
                      0.6,0.95, 
                      0.4,0.05,
                      0.1,0.35)
  }
  if (batch_size==6){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05)
  }
  if (batch_size==8){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6)
  }
  if (batch_size==10){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6,
                      
                      0.2,0.8,
                      0.8,0.1)
  }
  if (batch_size==5){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.4,0.05,
                      0.1,0.35)
  }
  if (batch_size==7){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.7,0.6)
    
  }
  if (batch_size==9){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6,
                      
                      0.2,0.8)
  }
  
  
  
  
  #this creates "number_ICs" lots of initial starting points
  #1st is with no noise
  #2nd onwards are with 2D multivariate normal noise mean (0,0) var=(sigma_noise^2 , 0
  #                                                                  0 , sigma_noise^2)
  noise_vector <- noise_function(batch_size,
                                 number_ICs,
                                 sigma_noise)
  
  starting_pts_without_ICs <- rep(starting_pts,number_ICs)
  starting_pts_with_ICs <- starting_pts_without_ICs+noise_vector
  
  
  
  
  #this ensures any noisy initial values below 0 are changed to 0.01
  #similarly above 1 changed to 0.99
  starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs<0,0.01)
  starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs>1,0.99)
  
  
  
  
  #now we find the optim value over all ICs
  #we save the parameters and the optimised value
  optim_par <- rep(0,number_ICs*batch_size*2)
  optim_value <- rep(0,number_ICs)
  for (j in 1:number_ICs){
    optim <- optim(method = "L-BFGS-B",
                   
                   par=starting_pts_with_ICs[(2*(j-1)*batch_size+1):(2*j*batch_size)], 
                   
                   fn=sigmoidal_optim_function,
                   
                   oldpts=batch0,
                   grid_length=grid_length,
                   theta_val2=theta_emulator,
                   sigma_val2=sigma_emulator,
                   E_f_val2=E_emulator,
                   nugget_val2=nugget_emulator,
                   true_runs_vec=batch0_evaluated,
                   
                   runtime_emulator_theta=run_time_emulator_theta,
                   runtime_emulator_sigma=run_time_emulator_sigma,
                   runtime_emulator_E=run_time_emulator_E,
                   runtime_emulator_nugget=run_time_emulator_nugget,
                   
                   ceiling=ceiling_value,
                   freeze=FALSE, #this step is just about doing the melt (gentle) sigmoidal
                   
                   lower=0,
                   upper=1)
    optim_par[((batch_size*2)*(j-1)+1):((batch_size*2)*j)] <- optim$par
    optim_value[j] <- optim$value
  }
  
  
  
  
  
  
  #this produces var graph showing locations of optimised points for every IC for every alpha
  #warning - this can produce hundreds of graphs!
  if (want_all_var_graphs==TRUE){
    for (z in 1:(number_ICs)){
      var_graph(xD_vec_original=batch0,
                xD_vec_new=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                grid_length=grid_length, 
                theta_val=theta_emulator, 
                sigma_val=sigma_emulator, 
                E_f_val=E_emulator,
                nugget_val=nugget_emulator)
    }
  }
  
  
  
  
  #we calculate v90 and Tc for each set of optimised locations
  #calculate V90 and Tc for each set of locations
  oldpts_matrix <- matrix(batch0,ncol=2,byrow=TRUE) 
  v90 <- rep(0,number_ICs)
  Tc <- rep(0,number_ICs)
  for (z in 1:(number_ICs)){ 
    newpts_matrix <- matrix(optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],ncol=2,byrow=TRUE)
    xD <- rbind(oldpts_matrix,newpts_matrix)
    D <- f(xD)
    x_grid <- seq(-0.001,1.001,len=grid_length)
    xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
    em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_emulator,sigma=sigma_emulator,E_f=E_emulator,nugget=nugget_emulator)
    v90[z] <- quantile(em_out[,"VarD_f(x)"],percentile)
    Tc[z] <- emulate_log_run_time_Tc(new_batch_vector=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                                     xD_matrix=oldpts_matrix,
                                     D_vector=batch0_evaluated,
                                     theta_val=run_time_emulator_theta,
                                     sigma_val=run_time_emulator_sigma,
                                     E_f_val=run_time_emulator_E,
                                     nugget_val=run_time_emulator_nugget) #take care this should use the log emulator characterisitics!!! May have to change in future
  }
  
  
  
  
  
  #best design (over the ICs) based on which gives the lowest optim value
  minimum_index <- which.min(optim_value)
  
  best_design_v90 <- v90[minimum_index]
  best_design_Tc <- Tc[minimum_index]
  best_design_locations <- matrix(optim_par[((batch_size*2)*(minimum_index-1)+1):((batch_size*2)*minimum_index)],ncol=2,byrow=TRUE)
  best_design_value <- optim_value[minimum_index]
  
  
  
  #plots the variance graph (showing where optimised locations are) only for the best design
  if (want_only_best_var_graph==TRUE){
    var_graph(xD_vec_original=batch0,
              xD_vec_new=c(t(best_design_locations)),
              grid_length=grid_length, 
              theta_val=theta_emulator, 
              sigma_val=sigma_emulator, 
              E_f_val=E_emulator,
              nugget_val=nugget_emulator)
    title("Best Design",adj=0,line=-26)
    
  }
  
  
  
  
  freeze_v90 <- NA
  freeze_Tc <- NA
  freeze_locations <- rep(NA,batch_size)
  freeze_value <- NA
  
  if (want_freeze==TRUE){
    freeze_optim <- optim(method = "L-BFGS-B",
                          
                          par=c(t(best_design_locations)), 
                          
                          fn=sigmoidal_optim_function,
                          
                          oldpts=batch0,
                          grid_length=grid_length,
                          theta_val2=theta_emulator,
                          sigma_val2=sigma_emulator,
                          E_f_val2=E_emulator,
                          nugget_val2=nugget_emulator,
                          true_runs_vec=batch0_evaluated,
                          
                          runtime_emulator_theta=run_time_emulator_theta,
                          runtime_emulator_sigma=run_time_emulator_sigma,
                          runtime_emulator_E=run_time_emulator_E,
                          runtime_emulator_nugget=run_time_emulator_nugget,
                          
                          ceiling=ceiling_value,
                          freeze=TRUE, #now we use a stricter sigmoidal!
                          
                          lower=0,
                          upper=1)
    freeze_locations <- freeze_optim$par
    freeze_value <- freeze_optim$value
    
    var_graph(xD_vec_original=batch0,
              xD_vec_new=freeze_locations,
              grid_length=grid_length, 
              theta_val=theta_emulator, 
              sigma_val=sigma_emulator, 
              E_f_val=E_emulator,
              nugget_val=nugget_emulator)
    title("Freeze Design",adj=0,line=-26)
    
    oldpts_matrix <- matrix(batch0,ncol=2,byrow=TRUE)
    freeze_locations_matrix <- matrix(freeze_locations,ncol=2,byrow=TRUE)
    xD <- rbind(oldpts_matrix,freeze_locations_matrix)
    D <- f(xD)
    x_grid <- seq(-0.001,1.001,len=grid_length)
    xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
    em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_emulator,sigma=sigma_emulator,E_f=E_emulator,nugget=nugget_emulator)
    freeze_v90 <- quantile(em_out[,"VarD_f(x)"],percentile)
    freeze_Tc <- emulate_log_run_time_Tc(new_batch_vector=freeze_locations,
                                         xD_matrix=oldpts_matrix,
                                         D_vector=batch0_evaluated,
                                         theta_val=run_time_emulator_theta,
                                         sigma_val=run_time_emulator_sigma,
                                         E_f_val=run_time_emulator_E,
                                         nugget_val=run_time_emulator_nugget)
  }
  
  
  
  #return some things we may want to look at
  return(list("All v90"=v90,
              "All Tc"=Tc,
              "All locations"=optim_par,
              "All values"=optim_value,
              
              "Best Design v90"=best_design_v90,
              "Best Design Tc"=best_design_Tc,
              "Best Design locations"=best_design_locations,
              "Best Design value"=best_design_value,
              
              "Freeze v90"=freeze_v90,
              "Freeze Tc"=freeze_Tc,
              "Freeze location"=freeze_locations,
              "Freeze value"=freeze_value,
              
              "Best design index"=minimum_index
  ))
}


## end of run before ceiling. ----------------------------------------------




# Figure 38 ---------------------------------------------------------------

## designs when theta=1 ----------------------------------------------------

### 4 -----------------------------------------------------------------------
#png(file="q_theta1_pts4_6000_%d.png",width=1400,height=1200,res=200)
q_theta1_pts4_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                    want_only_best_var_graph=TRUE,
                                                    want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                    
                                                    batch0=c(t(xD)),
                                                    batch0_evaluated=D, #evaluated run times
                                                    
                                                    grid_length=60,
                                                    
                                                    batch_size=4, #4, 6, 8 or 10
                                                    
                                                    theta_emulator=1,
                                                    sigma_emulator=sd(D),
                                                    E_emulator=mean(D),
                                                    nugget_emulator=0.001^2, #input as squared
                                                    
                                                    run_time_emulator_theta=1,
                                                    run_time_emulator_sigma=sd(D),
                                                    run_time_emulator_E=mean(D),
                                                    run_time_emulator_nugget=0.145064^2,
                                                    
                                                    percentile=0.9, #percentile of variance used in optim 
                                                    number_ICs=40, #how many different starting points we want
                                                    sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                    
                                                    ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()


### 5 -----------------------------------------------------------------------
#png(file="q_theta1_pts5_6000_%d.png",width=1400,height=1200,res=200)
q_theta1_pts5_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                    want_only_best_var_graph=TRUE,
                                                    want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                    
                                                    batch0=c(t(xD)),
                                                    batch0_evaluated=D, #evaluated run times
                                                    
                                                    grid_length=60,
                                                    
                                                    batch_size=5, #4, 6, 8 or 10
                                                    
                                                    theta_emulator=1,
                                                    sigma_emulator=sd(D),
                                                    E_emulator=mean(D),
                                                    nugget_emulator=0.001^2, #input as squared
                                                    
                                                    run_time_emulator_theta=1,
                                                    run_time_emulator_sigma=sd(D),
                                                    run_time_emulator_E=mean(D),
                                                    run_time_emulator_nugget=0.145064^2,
                                                    
                                                    percentile=0.9, #percentile of variance used in optim 
                                                    number_ICs=40, #how many different starting points we want
                                                    sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                    
                                                    ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()



### 6 -----------------------------------------------------------------------
#png(file="q_theta1_pts6_6000_%d.png",width=1400,height=1200,res=200)
q_theta1_pts6_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                    want_only_best_var_graph=TRUE,
                                                    want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                    
                                                    batch0=c(t(xD)),
                                                    batch0_evaluated=D, #evaluated run times
                                                    
                                                    grid_length=60,
                                                    
                                                    batch_size=6, #4, 6, 8 or 10
                                                    
                                                    theta_emulator=1,
                                                    sigma_emulator=sd(D),
                                                    E_emulator=mean(D),
                                                    nugget_emulator=0.001^2, #input as squared
                                                    
                                                    run_time_emulator_theta=1,
                                                    run_time_emulator_sigma=sd(D),
                                                    run_time_emulator_E=mean(D),
                                                    run_time_emulator_nugget=0.145064^2,
                                                    
                                                    percentile=0.9, #percentile of variance used in optim 
                                                    number_ICs=40, #how many different starting points we want
                                                    sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                    
                                                    ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()



### 7 -----------------------------------------------------------------------
#png(file="q_theta1_pts7_6000_%d.png",width=1400,height=1200,res=200)
q_theta1_pts7_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                    want_only_best_var_graph=TRUE,
                                                    want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                    
                                                    batch0=c(t(xD)),
                                                    batch0_evaluated=D, #evaluated run times
                                                    
                                                    grid_length=60,
                                                    
                                                    batch_size=7, #4, 6, 8 or 10
                                                    
                                                    theta_emulator=1,
                                                    sigma_emulator=sd(D),
                                                    E_emulator=mean(D),
                                                    nugget_emulator=0.001^2, #input as squared
                                                    
                                                    run_time_emulator_theta=1,
                                                    run_time_emulator_sigma=sd(D),
                                                    run_time_emulator_E=mean(D),
                                                    run_time_emulator_nugget=0.145064^2,
                                                    
                                                    percentile=0.9, #percentile of variance used in optim 
                                                    number_ICs=40, #how many different starting points we want
                                                    sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                    
                                                    ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()



### 8 -----------------------------------------------------------------------
#png(file="q_theta1_pts8_6000_%d.png",width=1400,height=1200,res=200)
q_theta1_pts8_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                    want_only_best_var_graph=TRUE,
                                                    want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                    
                                                    batch0=c(t(xD)),
                                                    batch0_evaluated=D, #evaluated run times
                                                    
                                                    grid_length=60,
                                                    
                                                    batch_size=8, #4, 6, 8 or 10
                                                    
                                                    theta_emulator=1,
                                                    sigma_emulator=sd(D),
                                                    E_emulator=mean(D),
                                                    nugget_emulator=0.001^2, #input as squared
                                                    
                                                    run_time_emulator_theta=1,
                                                    run_time_emulator_sigma=sd(D),
                                                    run_time_emulator_E=mean(D),
                                                    run_time_emulator_nugget=0.145064^2,
                                                    
                                                    percentile=0.9, #percentile of variance used in optim 
                                                    number_ICs=40, #how many different starting points we want
                                                    sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                    
                                                    ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()



### 9 -----------------------------------------------------------------------
#png(file="q_theta1_pts9_6000_%d.png",width=1400,height=1200,res=200)
q_theta1_pts9_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                    want_only_best_var_graph=TRUE,
                                                    want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                    
                                                    batch0=c(t(xD)),
                                                    batch0_evaluated=D, #evaluated run times
                                                    
                                                    grid_length=60,
                                                    
                                                    batch_size=9, #4, 6, 8 or 10
                                                    
                                                    theta_emulator=1,
                                                    sigma_emulator=sd(D),
                                                    E_emulator=mean(D),
                                                    nugget_emulator=0.001^2, #input as squared
                                                    
                                                    run_time_emulator_theta=1,
                                                    run_time_emulator_sigma=sd(D),
                                                    run_time_emulator_E=mean(D),
                                                    run_time_emulator_nugget=0.145064^2,
                                                    
                                                    percentile=0.9, #percentile of variance used in optim 
                                                    number_ICs=40, #how many different starting points we want
                                                    sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                    
                                                    ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()





# Figure 39 ---------------------------------------------------------------

## designs when theta=0.3 --------------------------------------------------


### 4 -----------------------------------------------------------------------

#png(file="q_thetapoint3_pts4_6000_%d.png",width=1400,height=1200,res=200)
q_thetapoint3_pts4_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                         want_only_best_var_graph=TRUE,
                                                         want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                         
                                                         batch0=c(t(xD)),
                                                         batch0_evaluated=D, #evaluated run times
                                                         
                                                         grid_length=60,
                                                         
                                                         batch_size=4, #4, 6, 8 or 10
                                                         
                                                         theta_emulator=0.3,
                                                         sigma_emulator=sd(D),
                                                         E_emulator=mean(D),
                                                         nugget_emulator=0.001^2, #input as squared
                                                         
                                                         run_time_emulator_theta=1,
                                                         run_time_emulator_sigma=sd(D),
                                                         run_time_emulator_E=mean(D),
                                                         run_time_emulator_nugget=0.145064^2,
                                                         
                                                         percentile=0.9, #percentile of variance used in optim 
                                                         number_ICs=40, #how many different starting points we want
                                                         sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                         
                                                         ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()

### 5 -----------------------------------------------------------------------

#png(file="q_thetapoint3_pts5_6000_%d.png",width=1400,height=1200,res=200)
q_thetapoint3_pts5_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                         want_only_best_var_graph=TRUE,
                                                         want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                         
                                                         batch0=c(t(xD)),
                                                         batch0_evaluated=D, #evaluated run times
                                                         
                                                         grid_length=60,
                                                         
                                                         batch_size=5, #4, 6, 8 or 10
                                                         
                                                         theta_emulator=0.3,
                                                         sigma_emulator=sd(D),
                                                         E_emulator=mean(D),
                                                         nugget_emulator=0.001^2, #input as squared
                                                         
                                                         run_time_emulator_theta=1,
                                                         run_time_emulator_sigma=sd(D),
                                                         run_time_emulator_E=mean(D),
                                                         run_time_emulator_nugget=0.145064^2,
                                                         
                                                         percentile=0.9, #percentile of variance used in optim 
                                                         number_ICs=40, #how many different starting points we want
                                                         sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                         
                                                         ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()


### 6 -----------------------------------------------------------------------

#png(file="q_thetapoint3_pts6_6000_%d.png",width=1400,height=1200,res=200)
q_thetapoint3_pts6_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                         want_only_best_var_graph=TRUE,
                                                         want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                         
                                                         batch0=c(t(xD)),
                                                         batch0_evaluated=D, #evaluated run times
                                                         
                                                         grid_length=60,
                                                         
                                                         batch_size=6, #4, 6, 8 or 10
                                                         
                                                         theta_emulator=0.3,
                                                         sigma_emulator=sd(D),
                                                         E_emulator=mean(D),
                                                         nugget_emulator=0.001^2, #input as squared
                                                         
                                                         run_time_emulator_theta=1,
                                                         run_time_emulator_sigma=sd(D),
                                                         run_time_emulator_E=mean(D),
                                                         run_time_emulator_nugget=0.145064^2,
                                                         
                                                         percentile=0.9, #percentile of variance used in optim 
                                                         number_ICs=40, #how many different starting points we want
                                                         sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                         
                                                         ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()


### 7 -----------------------------------------------------------------------

#png(file="q_thetapoint3_pts7_6000_%d.png",width=1400,height=1200,res=200)
q_thetapoint3_pts7_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                         want_only_best_var_graph=TRUE,
                                                         want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                         
                                                         batch0=c(t(xD)),
                                                         batch0_evaluated=D, #evaluated run times
                                                         
                                                         grid_length=60,
                                                         
                                                         batch_size=7, #4, 6, 8 or 10
                                                         
                                                         theta_emulator=0.3,
                                                         sigma_emulator=sd(D),
                                                         E_emulator=mean(D),
                                                         nugget_emulator=0.001^2, #input as squared
                                                         
                                                         run_time_emulator_theta=1,
                                                         run_time_emulator_sigma=sd(D),
                                                         run_time_emulator_E=mean(D),
                                                         run_time_emulator_nugget=0.145064^2,
                                                         
                                                         percentile=0.9, #percentile of variance used in optim 
                                                         number_ICs=40, #how many different starting points we want
                                                         sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                         
                                                         ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()


### 8 -----------------------------------------------------------------------

#png(file="q_thetapoint3_pts8_6000_%d.png",width=1400,height=1200,res=200)
q_thetapoint3_pts8_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                         want_only_best_var_graph=TRUE,
                                                         want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                         
                                                         batch0=c(t(xD)),
                                                         batch0_evaluated=D, #evaluated run times
                                                         
                                                         grid_length=60,
                                                         
                                                         batch_size=8, #4, 6, 8 or 10
                                                         
                                                         theta_emulator=0.3,
                                                         sigma_emulator=sd(D),
                                                         E_emulator=mean(D),
                                                         nugget_emulator=0.001^2, #input as squared
                                                         
                                                         run_time_emulator_theta=1,
                                                         run_time_emulator_sigma=sd(D),
                                                         run_time_emulator_E=mean(D),
                                                         run_time_emulator_nugget=0.145064^2,
                                                         
                                                         percentile=0.9, #percentile of variance used in optim 
                                                         number_ICs=40, #how many different starting points we want
                                                         sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                         
                                                         ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()


### 9 -----------------------------------------------------------------------

#png(file="q_thetapoint3_pts9_6000_%d.png",width=1400,height=1200,res=200)
q_thetapoint3_pts9_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                         want_only_best_var_graph=TRUE,
                                                         want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                         
                                                         batch0=c(t(xD)),
                                                         batch0_evaluated=D, #evaluated run times
                                                         
                                                         grid_length=60,
                                                         
                                                         batch_size=9, #4, 6, 8 or 10
                                                         
                                                         theta_emulator=0.3,
                                                         sigma_emulator=sd(D),
                                                         E_emulator=mean(D),
                                                         nugget_emulator=0.001^2, #input as squared
                                                         
                                                         run_time_emulator_theta=1,
                                                         run_time_emulator_sigma=sd(D),
                                                         run_time_emulator_E=mean(D),
                                                         run_time_emulator_nugget=0.145064^2,
                                                         
                                                         percentile=0.9, #percentile of variance used in optim 
                                                         number_ICs=40, #how many different starting points we want
                                                         sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                         
                                                         ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()






# Figure 40 ---------------------------------------------------------------

## run this! ---------------------------------------------------------------
rm(list = ls())


library(viridisLite)
library(pdist)

nugget_argument_efficient_nugget_all_graphs <- function(xD_matrix, 
                                                        D_vector,
                                                        grid_length=50, 
                                                        theta_val=0.5, 
                                                        sigma_val=0.5, 
                                                        E_f_val=0,
                                                        
                                                        #new part:
                                                        nugget_val=10^(-6)){
  
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
  #new part (adding nugget as argument):
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)  
  
  
  #extracting the emulator expectation and variance
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  
  #graphing
  emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")
}





nugget_argument_efficient_nugget_simple_BL_emulator <- function(xP,             # the set of emulator prediction points
                                                                xD,             # the run input locations xD
                                                                D,              # the run outputs D = (f(x^1),...,f(x^n))
                                                                theta = 1,      # the correlation lengths (can be a vector)
                                                                sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                                                E_f = 0,        # prior expectation of f: E(f(x)) = 0 
                                                                using_pdist = 1,# if you have installed pdist package
                                                                
                                                                #new part:
                                                                nugget = 10^(-6) #must input nugget 0.1^2 
                                                                #squared because working with sigma squared
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
  
  
  
  
  
  
  #deterministic: manage numerical error
  #diagonal <- diag(10^-6,n,n) #1 millionth of what we are doing but well above numerical error (nice safe amount)
  #stochastic: manage stochasticity
  #diagonal <- diag(0.1^2,n,n) #each green pt quite a bit more uncertain
  
  #new part:
  diagonal <- diag(nugget,n,n)
  Var_D <- Var_D + sigma^2*diagonal
  
  #note that nugget 0.1^2 or 0.2^2 closer to 10% of mean of D (which is roughly what we want)
  #could put in directly but will have to change by hand every time change emulator parameters
  
  
  
  
  
  
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



dataframe_sim_run_times <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/simulation_overview.txt")
matrix_sim_run_times <- as.matrix(dataframe_sim_run_times)
colnames(matrix_sim_run_times) <- c("Emulator Coordinate x", 
                                    "Emulator Coordinate y", 
                                    "Runtime [CPU-Hours]", 
                                    "Batch Number", 
                                    "fmin",	
                                    "fmax",	
                                    "rho0 [cm^-3]",
                                    "T_AGN [K]")
x <- matrix_sim_run_times[,"Emulator Coordinate x"]
y <- matrix_sim_run_times[,"Emulator Coordinate y"]
plot(x,y)

xD <- matrix(data=c(x,y),ncol=2,byrow = FALSE)
D <- log10(matrix_sim_run_times[,"Runtime [CPU-Hours]"])

emulate_log_run_time_Tc <- function(new_batch_vector,
                                    xD_matrix,
                                    D_vector,
                                    
                                    theta_val,
                                    sigma_val,
                                    E_f_val,
                                    
                                    nugget_val=10^(-6))
  
{
  xD <- xD_matrix
  D <- D_vector
  new_batch_matrix <- matrix(data=new_batch_vector,ncol=2,byrow = TRUE)
  xP <- new_batch_matrix
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)
  E_D_fx_mat <- em_out[,"ExpD_f(x)"] 
  normal_times <- 10^(E_D_fx_mat)
  return(sum(normal_times))
}


var_graph <- function(xD_vec_original,
                      xD_vec_new,
                      grid_length=50, 
                      theta_val=0.5, 
                      sigma_val=0.5, 
                      E_f_val=0,
                      
                      nugget_val=10^(-6))
{
  par(mar=c(4,4,3,3))
  var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
  
  xD_vec <- c(xD_vec_original,xD_vec_new)
  xD <- matrix(xD_vec,ncol=2,byrow=TRUE)
  D <- f(xD)
  
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)   
  
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]",
                 #new part:
                 xD_col=rep(c("purple","pink"),c(length(xD_vec_original)/2,length(xD_vec_new)/2)) #/2 because x and y coords in vector so number of points is length/2
  )
  
}

#create random function (won't affect variance but function needs it)
f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))





library(MASS)

noise_function <- function(batch_size,
                           number_ICs,
                           sigma_noise){
  var_matrix <- matrix(c(sigma_noise^2,0,
                         0,sigma_noise^2),
                       nrow=2,
                       ncol=2,
                       byrow=TRUE
  )
  noise2D_matrix <- mvrnorm(n=batch_size*(number_ICs-1),
                            mu=rep(0,2),
                            Sigma=var_matrix)
  noise2D_vector <- c(rep(0,batch_size*2),c(t(noise2D_matrix)))
  return(noise2D_vector)
}


# This function had a slight flaw
# Previously we calculated 0.9 (90th percentile) every time
# Here we need 0.95 therefore we change 0.9 to percentile3 
# in sigmoidal_optim_function below
# this flaw does not affect any other plots (as we use 0.9 everywhere except here)
# useful to note

sigmoidal_optim_function <- function(newpts, #will optimise first argument
                                     oldpts,
                                     grid_length=50,
                                     theta_val2=0.5,
                                     sigma_val2=0.5,
                                     E_f_val2=0,
                                     nugget_val2=10^(-6),
                                     true_runs_vec,
                                     
                                     ceiling,
                                     freeze,
                                     
                                     runtime_emulator_theta,
                                     runtime_emulator_sigma,
                                     runtime_emulator_E,
                                     runtime_emulator_nugget,
                                     
                                     #NEW!
                                     percentile_val3
){
  oldpts_matrix <- matrix(oldpts,ncol=2,byrow=TRUE) 
  newpts_matrix <- matrix(newpts,ncol=2,byrow=TRUE)
  xD <- rbind(oldpts_matrix,newpts_matrix)
  D <- f(xD)
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val2,sigma=sigma_val2,E_f=E_f_val2,nugget=nugget_val2)
  #CHANGED THIS LINE!!!
  v90 <- quantile(em_out[,"VarD_f(x)"],percentile_val3)
  Tc <- emulate_log_run_time_Tc(new_batch_vector=newpts,
                                xD_matrix=oldpts_matrix,
                                D_vector=true_runs_vec,
                                theta_val=runtime_emulator_theta,
                                sigma_val=runtime_emulator_sigma,
                                E_f_val=runtime_emulator_E,
                                nugget_val=runtime_emulator_nugget)
  
  if (freeze==FALSE){
    a <- 5 #steepness
    b <- 3 #height
    c <- ceiling #log10(ceiling) x value halfway up the sigmoidal
    d <- 0.5 #steepness of linear term added on 
  }
  
  if (freeze==TRUE){
    a <- 10 #steepness - THIS IS THE ONLY DIFFERENCE
    b <- 3 #height
    c <- ceiling #log10(ceiling) x value halfway up the sigmoidal
    d <- 0.5 #steepness of linear term added on
  }
  
  x <- Tc
  
  if (Tc<=ceiling){
    sigmoidal <- b/(1+exp(-a*(x-c)))
  }
  if (Tc>ceiling){
    sigmoidal <- b/(1+exp(-a*(x-c))) + d*x - d*c
  }
  
  final <- (1/2)*log10(v90) + sigmoidal
  return(final)
}

# #checking sigmoidal part works ok
# sig <- function(log10Tc,
#                 log10ceiling
#                 ){
#   a <- 10 #5 or 10
#   b <- 3
#   c <- log10ceiling
#   d <- 0.5
#   x <- log10Tc
#   if (log10Tc<=log10ceiling){
#     sigmoidal <- b/(1+exp(-a*(x-c)))
#   }
#   if (log10Tc>log10ceiling){
#     sigmoidal <- b/(1+exp(-a*(x-c))) + d*x - d*c
#   }
#   return(sigmoidal)
# }
# x <- seq(0,10,by=0.01)
# which(x==4)
# plot(x,sapply(x,FUN=sig,log10ceiling=7),ylim=c(0,6))
# y <- sapply(x,FUN=sig,log10ceiling=7)
# y[601]




#ALSO HAVE TO CHANGE THIS SLIGHTLY NOW


#Function using sigmoidal criteria
everything_function_sigmoidal <- function(want_all_var_graphs,
                                          want_only_best_var_graph,
                                          want_freeze, #initialise design at melt design found now with harsher sigmoidal
                                          
                                          batch0,
                                          batch0_evaluated, #evaluated run times
                                          
                                          grid_length,
                                          
                                          batch_size, #4, 6, 8 or 10 now 5,7 also allowed
                                          
                                          theta_emulator,
                                          sigma_emulator,
                                          E_emulator,
                                          nugget_emulator, #input as squared
                                          
                                          run_time_emulator_theta,
                                          run_time_emulator_sigma,
                                          run_time_emulator_E,
                                          run_time_emulator_nugget,
                                          
                                          percentile, #percentile of variance used in optim 
                                          number_ICs, #how many different starting points we want
                                          sigma_noise, #N(0,sigma_noise^2) for random noise
                                          
                                          ceiling_value #NOT LOGGED max total run time cannot exceed
){
  
  #create random function (won't affect variance but function needs it)
  f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))
  
  
  
  #determine initial starting points 
  #dependent on batch_size (number of new points) we want
  if (batch_size==4){
    starting_pts <- c(0.9,0.85,
                      0.6,0.95, 
                      0.4,0.05,
                      0.1,0.35)
  }
  if (batch_size==6){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05)
  }
  if (batch_size==8){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6)
  }
  if (batch_size==10){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6,
                      
                      0.2,0.8,
                      0.8,0.1)
  }
  if (batch_size==5){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.4,0.05,
                      0.1,0.35)
  }
  if (batch_size==7){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.7,0.6)
    
  }
  if (batch_size==9){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6,
                      
                      0.2,0.8)
  }
  
  
  
  
  #this creates "number_ICs" lots of initial starting points
  #1st is with no noise
  #2nd onwards are with 2D multivariate normal noise mean (0,0) var=(sigma_noise^2 , 0
  #                                                                  0 , sigma_noise^2)
  noise_vector <- noise_function(batch_size,
                                 number_ICs,
                                 sigma_noise)
  
  starting_pts_without_ICs <- rep(starting_pts,number_ICs)
  starting_pts_with_ICs <- starting_pts_without_ICs+noise_vector
  
  
  
  
  #this ensures any noisy initial values below 0 are changed to 0.01
  #similarly above 1 changed to 0.99
  starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs<0,0.01)
  starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs>1,0.99)
  
  
  
  
  #now we find the optim value over all ICs
  #we save the parameters and the optimised value
  optim_par <- rep(0,number_ICs*batch_size*2)
  optim_value <- rep(0,number_ICs)
  for (j in 1:number_ICs){
    optim <- optim(method = "L-BFGS-B",
                   
                   par=starting_pts_with_ICs[(2*(j-1)*batch_size+1):(2*j*batch_size)], 
                   
                   fn=sigmoidal_optim_function,
                   
                   oldpts=batch0,
                   grid_length=grid_length,
                   theta_val2=theta_emulator,
                   sigma_val2=sigma_emulator,
                   E_f_val2=E_emulator,
                   nugget_val2=nugget_emulator,
                   true_runs_vec=batch0_evaluated,
                   
                   runtime_emulator_theta=run_time_emulator_theta,
                   runtime_emulator_sigma=run_time_emulator_sigma,
                   runtime_emulator_E=run_time_emulator_E,
                   runtime_emulator_nugget=run_time_emulator_nugget,
                   
                   ceiling=ceiling_value,
                   freeze=FALSE, #this step is just about doing the melt (gentle) sigmoidal
                   
                   percentile_val3 = percentile,
                   
                   lower=0,
                   upper=1)
    optim_par[((batch_size*2)*(j-1)+1):((batch_size*2)*j)] <- optim$par
    optim_value[j] <- optim$value
  }
  
  
  
  
  
  
  #this produces var graph showing locations of optimised points for every IC for every alpha
  #warning - this can produce hundreds of graphs!
  if (want_all_var_graphs==TRUE){
    for (z in 1:(number_ICs)){
      var_graph(xD_vec_original=batch0,
                xD_vec_new=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                grid_length=grid_length, 
                theta_val=theta_emulator, 
                sigma_val=sigma_emulator, 
                E_f_val=E_emulator,
                nugget_val=nugget_emulator)
    }
  }
  
  
  
  
  #we calculate v90 and Tc for each set of optimised locations
  #calculate V90 and Tc for each set of locations
  oldpts_matrix <- matrix(batch0,ncol=2,byrow=TRUE) 
  v90 <- rep(0,number_ICs)
  Tc <- rep(0,number_ICs)
  for (z in 1:(number_ICs)){ 
    newpts_matrix <- matrix(optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],ncol=2,byrow=TRUE)
    xD <- rbind(oldpts_matrix,newpts_matrix)
    D <- f(xD)
    x_grid <- seq(-0.001,1.001,len=grid_length)
    xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
    em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_emulator,sigma=sigma_emulator,E_f=E_emulator,nugget=nugget_emulator)
    v90[z] <- quantile(em_out[,"VarD_f(x)"],percentile)
    Tc[z] <- emulate_log_run_time_Tc(new_batch_vector=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                                     xD_matrix=oldpts_matrix,
                                     D_vector=batch0_evaluated,
                                     theta_val=run_time_emulator_theta,
                                     sigma_val=run_time_emulator_sigma,
                                     E_f_val=run_time_emulator_E,
                                     nugget_val=run_time_emulator_nugget) #take care this should use the log emulator characterisitics!!! May have to change in future
  }
  
  
  
  
  
  #best design (over the ICs) based on which gives the lowest optim value
  minimum_index <- which.min(optim_value)
  
  best_design_v90 <- v90[minimum_index]
  best_design_Tc <- Tc[minimum_index]
  best_design_locations <- matrix(optim_par[((batch_size*2)*(minimum_index-1)+1):((batch_size*2)*minimum_index)],ncol=2,byrow=TRUE)
  best_design_value <- optim_value[minimum_index]
  
  
  
  #plots the variance graph (showing where optimised locations are) only for the best design
  if (want_only_best_var_graph==TRUE){
    var_graph(xD_vec_original=batch0,
              xD_vec_new=c(t(best_design_locations)),
              grid_length=grid_length, 
              theta_val=theta_emulator, 
              sigma_val=sigma_emulator, 
              E_f_val=E_emulator,
              nugget_val=nugget_emulator)
    title("Best Design",adj=0,line=-26)
    
  }
  
  
  
  
  freeze_v90 <- NA
  freeze_Tc <- NA
  freeze_locations <- rep(NA,batch_size)
  freeze_value <- NA
  
  if (want_freeze==TRUE){
    freeze_optim <- optim(method = "L-BFGS-B",
                          
                          par=c(t(best_design_locations)), 
                          
                          fn=sigmoidal_optim_function,
                          
                          oldpts=batch0,
                          grid_length=grid_length,
                          theta_val2=theta_emulator,
                          sigma_val2=sigma_emulator,
                          E_f_val2=E_emulator,
                          nugget_val2=nugget_emulator,
                          true_runs_vec=batch0_evaluated,
                          
                          runtime_emulator_theta=run_time_emulator_theta,
                          runtime_emulator_sigma=run_time_emulator_sigma,
                          runtime_emulator_E=run_time_emulator_E,
                          runtime_emulator_nugget=run_time_emulator_nugget,
                          
                          ceiling=ceiling_value,
                          freeze=TRUE, #now we use a stricter sigmoidal!
                          
                          lower=0,
                          upper=1)
    freeze_locations <- freeze_optim$par
    freeze_value <- freeze_optim$value
    
    var_graph(xD_vec_original=batch0,
              xD_vec_new=freeze_locations,
              grid_length=grid_length, 
              theta_val=theta_emulator, 
              sigma_val=sigma_emulator, 
              E_f_val=E_emulator,
              nugget_val=nugget_emulator)
    title("Freeze Design",adj=0,line=-26)
    
    oldpts_matrix <- matrix(batch0,ncol=2,byrow=TRUE)
    freeze_locations_matrix <- matrix(freeze_locations,ncol=2,byrow=TRUE)
    xD <- rbind(oldpts_matrix,freeze_locations_matrix)
    D <- f(xD)
    x_grid <- seq(-0.001,1.001,len=grid_length)
    xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
    em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_emulator,sigma=sigma_emulator,E_f=E_emulator,nugget=nugget_emulator)
    freeze_v90 <- quantile(em_out[,"VarD_f(x)"],percentile)
    freeze_Tc <- emulate_log_run_time_Tc(new_batch_vector=freeze_locations,
                                         xD_matrix=oldpts_matrix,
                                         D_vector=batch0_evaluated,
                                         theta_val=run_time_emulator_theta,
                                         sigma_val=run_time_emulator_sigma,
                                         E_f_val=run_time_emulator_E,
                                         nugget_val=run_time_emulator_nugget)
  }
  
  
  
  #return some things we may want to look at
  return(list("All v90"=v90,
              "All Tc"=Tc,
              "All locations"=optim_par,
              "All values"=optim_value,
              
              "Best Design v90"=best_design_v90,
              "Best Design Tc"=best_design_Tc,
              "Best Design locations"=best_design_locations,
              "Best Design value"=best_design_value,
              
              "Freeze v90"=freeze_v90,
              "Freeze Tc"=freeze_Tc,
              "Freeze location"=freeze_locations,
              "Freeze value"=freeze_value,
              
              "Best design index"=minimum_index
  ))
}




## now I will run V95 criteria! -------------------------------------

### 4 -----------------------------------------------------------------------

# png(file="q_v95_pts4_6000_%d.png",width=1400,height=1200,res=200)
# q_v95_pts4_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
#                                                  want_only_best_var_graph=TRUE,
#                                                  want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
# 
#                                                  batch0=c(t(xD)),
#                                                  batch0_evaluated=D, #evaluated run times
# 
#                                                  grid_length=60,
# 
#                                                  batch_size=4, #4, 6, 8 or 10
# 
#                                                  theta_emulator=1/2,
#                                                  sigma_emulator=sd(D),
#                                                  E_emulator=mean(D),
#                                                  nugget_emulator=0.001^2, #input as squared
# 
#                                                  run_time_emulator_theta=1,
#                                                  run_time_emulator_sigma=sd(D),
#                                                  run_time_emulator_E=mean(D),
#                                                  run_time_emulator_nugget=0.145064^2,
# 
#                                                  percentile=0.95, #percentile of variance used in optim
#                                                  number_ICs=40, #how many different starting points we want
#                                                  sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
# 
#                                                  ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
# )
# dev.off()

### 5 -----------------------------------------------------------------------

# png(file="q_v95_pts5_6000_%d.png",width=1400,height=1200,res=200)
# q_v95_pts5_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
#                                                  want_only_best_var_graph=TRUE,
#                                                  want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
# 
#                                                  batch0=c(t(xD)),
#                                                  batch0_evaluated=D, #evaluated run times
# 
#                                                  grid_length=60,
# 
#                                                  batch_size=5, #4, 6, 8 or 10
# 
#                                                  theta_emulator=1/2,
#                                                  sigma_emulator=sd(D),
#                                                  E_emulator=mean(D),
#                                                  nugget_emulator=0.001^2, #input as squared
# 
#                                                  run_time_emulator_theta=1,
#                                                  run_time_emulator_sigma=sd(D),
#                                                  run_time_emulator_E=mean(D),
#                                                  run_time_emulator_nugget=0.145064^2,
# 
#                                                  percentile=0.95, #percentile of variance used in optim
#                                                  number_ICs=40, #how many different starting points we want
#                                                  sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
# 
#                                                  ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
# )
# dev.off()
### 6 -----------------------------------------------------------------------

# png(file="q_v95_pts6_6000_%d.png",width=1400,height=1200,res=200)
# q_v95_pts6_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
#                                                  want_only_best_var_graph=TRUE,
#                                                  want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
# 
#                                                  batch0=c(t(xD)),
#                                                  batch0_evaluated=D, #evaluated run times
# 
#                                                  grid_length=60,
# 
#                                                  batch_size=6, #4, 6, 8 or 10
# 
#                                                  theta_emulator=1/2,
#                                                  sigma_emulator=sd(D),
#                                                  E_emulator=mean(D),
#                                                  nugget_emulator=0.001^2, #input as squared
# 
#                                                  run_time_emulator_theta=1,
#                                                  run_time_emulator_sigma=sd(D),
#                                                  run_time_emulator_E=mean(D),
#                                                  run_time_emulator_nugget=0.145064^2,
# 
#                                                  percentile=0.95, #percentile of variance used in optim
#                                                  number_ICs=40, #how many different starting points we want
#                                                  sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
# 
#                                                  ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
# )
# dev.off()
### 7 -----------------------------------------------------------------------

# png(file="q_v95_pts7_6000_%d.png",width=1400,height=1200,res=200)
# q_v95_pts7_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
#                                                  want_only_best_var_graph=TRUE,
#                                                  want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
# 
#                                                  batch0=c(t(xD)),
#                                                  batch0_evaluated=D, #evaluated run times
# 
#                                                  grid_length=60,
# 
#                                                  batch_size=7, #4, 6, 8 or 10
# 
#                                                  theta_emulator=1/2,
#                                                  sigma_emulator=sd(D),
#                                                  E_emulator=mean(D),
#                                                  nugget_emulator=0.001^2, #input as squared
# 
#                                                  run_time_emulator_theta=1,
#                                                  run_time_emulator_sigma=sd(D),
#                                                  run_time_emulator_E=mean(D),
#                                                  run_time_emulator_nugget=0.145064^2,
# 
#                                                  percentile=0.95, #percentile of variance used in optim
#                                                  number_ICs=40, #how many different starting points we want
#                                                  sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
# 
#                                                  ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
# )
# dev.off()
### 8 -----------------------------------------------------------------------

# png(file="q_v95_pts8_6000_%d.png",width=1400,height=1200,res=200)
# q_v95_pts8_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
#                                                  want_only_best_var_graph=TRUE,
#                                                  want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
# 
#                                                  batch0=c(t(xD)),
#                                                  batch0_evaluated=D, #evaluated run times
# 
#                                                  grid_length=60,
# 
#                                                  batch_size=8, #4, 6, 8 or 10
# 
#                                                  theta_emulator=1/2,
#                                                  sigma_emulator=sd(D),
#                                                  E_emulator=mean(D),
#                                                  nugget_emulator=0.001^2, #input as squared
# 
#                                                  run_time_emulator_theta=1,
#                                                  run_time_emulator_sigma=sd(D),
#                                                  run_time_emulator_E=mean(D),
#                                                  run_time_emulator_nugget=0.145064^2,
# 
#                                                  percentile=0.95, #percentile of variance used in optim
#                                                  number_ICs=40, #how many different starting points we want
#                                                  sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
# 
#                                                  ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
# )
# dev.off()
### 9 -----------------------------------------------------------------------

# png(file="q_v95_pts9_6000_%d.png",width=1400,height=1200,res=200)
# q_v95_pts9_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
#                                                  want_only_best_var_graph=TRUE,
#                                                  want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
# 
#                                                  batch0=c(t(xD)),
#                                                  batch0_evaluated=D, #evaluated run times
# 
#                                                  grid_length=60,
# 
#                                                  batch_size=9, #4, 6, 8 or 10
# 
#                                                  theta_emulator=1/2,
#                                                  sigma_emulator=sd(D),
#                                                  E_emulator=mean(D),
#                                                  nugget_emulator=0.001^2, #input as squared
# 
#                                                  run_time_emulator_theta=1,
#                                                  run_time_emulator_sigma=sd(D),
#                                                  run_time_emulator_E=mean(D),
#                                                  run_time_emulator_nugget=0.145064^2,
# 
#                                                  percentile=0.95, #percentile of variance used in optim
#                                                  number_ICs=40, #how many different starting points we want
#                                                  sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
# 
#                                                  ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
# )
# dev.off()






# Figure 41 ---------------------------------------------------------------

## run before alpha --------------------------------------------------------




rm(list = ls())

library(viridisLite)
library(pdist)

nugget_argument_efficient_nugget_all_graphs <- function(xD_matrix, 
                                                        D_vector,
                                                        grid_length=50, 
                                                        theta_val=0.5, 
                                                        sigma_val=0.5, 
                                                        E_f_val=0,
                                                        
                                                        #new part:
                                                        nugget_val=10^(-6)){
  
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
  #new part (adding nugget as argument):
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)  
  
  
  #extracting the emulator expectation and variance
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  
  #graphing
  emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")
}





nugget_argument_efficient_nugget_simple_BL_emulator <- function(xP,             # the set of emulator prediction points
                                                                xD,             # the run input locations xD
                                                                D,              # the run outputs D = (f(x^1),...,f(x^n))
                                                                theta = 1,      # the correlation lengths (can be a vector)
                                                                sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                                                E_f = 0,        # prior expectation of f: E(f(x)) = 0 
                                                                using_pdist = 1,# if you have installed pdist package
                                                                
                                                                #new part:
                                                                nugget = 10^(-6) #must input nugget 0.1^2 
                                                                #squared because working with sigma squared
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
  
  
  
  
  
  
  #deterministic: manage numerical error
  #diagonal <- diag(10^-6,n,n) #1 millionth of what we are doing but well above numerical error (nice safe amount)
  #stochastic: manage stochasticity
  #diagonal <- diag(0.1^2,n,n) #each green pt quite a bit more uncertain
  
  #new part:
  diagonal <- diag(nugget,n,n)
  Var_D <- Var_D + sigma^2*diagonal
  
  #note that nugget 0.1^2 or 0.2^2 closer to 10% of mean of D (which is roughly what we want)
  #could put in directly but will have to change by hand every time change emulator parameters
  
  
  
  
  
  
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


dataframe_sim_run_times <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/simulation_overview.txt")
matrix_sim_run_times <- as.matrix(dataframe_sim_run_times)
colnames(matrix_sim_run_times) <- c("Emulator Coordinate x", 
                                    "Emulator Coordinate y", 
                                    "Runtime [CPU-Hours]", 
                                    "Batch Number", 
                                    "fmin",	
                                    "fmax",	
                                    "rho0 [cm^-3]",
                                    "T_AGN [K]")
x <- matrix_sim_run_times[,"Emulator Coordinate x"]
y <- matrix_sim_run_times[,"Emulator Coordinate y"]

xD <- matrix(data=c(x,y),ncol=2,byrow = FALSE)
D <- log10(matrix_sim_run_times[,"Runtime [CPU-Hours]"])

emulate_log_run_time_Tc <- function(new_batch_vector,
                                    xD_matrix,
                                    D_vector,
                                    
                                    theta_val,
                                    sigma_val,
                                    E_f_val,
                                    
                                    nugget_val=10^(-6))
  
{
  xD <- xD_matrix
  D <- D_vector
  new_batch_matrix <- matrix(data=new_batch_vector,ncol=2,byrow = TRUE)
  xP <- new_batch_matrix
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)
  E_D_fx_mat <- em_out[,"ExpD_f(x)"] 
  normal_times <- 10^(E_D_fx_mat)
  return(sum(normal_times))
}

var_graph <- function(xD_vec_original,
                      xD_vec_new,
                      grid_length=50, 
                      theta_val=0.5, 
                      sigma_val=0.5, 
                      E_f_val=0,
                      
                      nugget_val=10^(-6))
{
  par(mar=c(4,4,3,3))
  var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
  
  xD_vec <- c(xD_vec_original,xD_vec_new)
  xD <- matrix(xD_vec,ncol=2,byrow=TRUE)
  D <- f(xD)
  
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)   
  
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]",
                 #new part:
                 xD_col=rep(c("purple","pink"),c(length(xD_vec_original)/2,length(xD_vec_new)/2)) #/2 because x and y coords in vector so number of points is length/2
  )
  
}

#create random function (won't affect variance but function needs it)
f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))


library(MASS)

noise_function <- function(batch_size,
                           number_ICs,
                           sigma_noise){
  var_matrix <- matrix(c(sigma_noise^2,0,
                         0,sigma_noise^2),
                       nrow=2,
                       ncol=2,
                       byrow=TRUE
  )
  noise2D_matrix <- mvrnorm(n=batch_size*(number_ICs-1),
                            mu=rep(0,2),
                            Sigma=var_matrix)
  noise2D_vector <- c(rep(0,batch_size*2),c(t(noise2D_matrix)))
  return(noise2D_vector)
}

new_optim_fixed_runs <- function(newpts, #will optimise first argument
                                 oldpts,
                                 grid_length=50,
                                 theta_val2=0.5,
                                 sigma_val2=0.5,
                                 E_f_val2=0,
                                 nugget_val2=10^(-6),
                                 true_runs_vec,
                                 alpha=1,
                                 
                                 runtime_emulator_theta,
                                 runtime_emulator_sigma,
                                 runtime_emulator_E,
                                 runtime_emulator_nugget
){ 
  oldpts_matrix <- matrix(oldpts,ncol=2,byrow=TRUE) 
  newpts_matrix <- matrix(newpts,ncol=2,byrow=TRUE)
  xD <- rbind(oldpts_matrix,newpts_matrix)
  D <- f(xD)
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val2,sigma=sigma_val2,E_f=E_f_val2,nugget=nugget_val2)
  v90 <- quantile(em_out[,"VarD_f(x)"],0.9)
  Tc <- emulate_log_run_time_Tc(new_batch_vector=newpts,
                                xD_matrix=oldpts_matrix,
                                D_vector=true_runs_vec,
                                theta_val=runtime_emulator_theta,
                                sigma_val=runtime_emulator_sigma,
                                E_f_val=runtime_emulator_E,
                                nugget_val=runtime_emulator_nugget)
  final <- (1/2)*log10(v90) + alpha*log10(Tc)
  return(final)
  
}

new_var_grid_graph <- function(want_var_graph,
                               want_grid,
                               want_graph,
                               want_best_var_graphs,
                               
                               batch0,
                               batch0_evaluated, #evaluated run times
                               
                               grid_length,
                               
                               batch_size, #4, 6, 8 or 10
                               
                               theta_emulator,
                               sigma_emulator,
                               E_emulator,
                               nugget_emulator, #input as squared
                               
                               run_time_emulator_theta,
                               run_time_emulator_sigma,
                               run_time_emulator_E,
                               run_time_emulator_nugget,
                               
                               alphas, #vector of alphas we want
                               percentile, #percentile of variance used in optim 
                               number_ICs, #how many different starting points we want
                               sigma_noise #N(0,sigma_noise^2) for random noise
){
  
  
  
  #determine initial starting points 
  #dependent on batch_size (number of new points) we want
  if (batch_size==4){
    starting_pts <- c(0.9,0.85,
                      0.6,0.95, 
                      0.4,0.05,
                      0.1,0.35)
  }
  if (batch_size==6){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05)
  }
  if (batch_size==8){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6)
  }
  if (batch_size==10){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6,
                      
                      0.2,0.8,
                      0.8,0.1)
  }
  
  
  
  
  #this creates "number_ICs" lots of initial starting points
  #1st is with no noise
  #2nd onwards are with 2D multivariate normal noise mean (0,0) var=(sigma_noise^2 , 0
  #                                                                  0 , sigma_noise^2)
  noise_vector <- noise_function(batch_size,
                                 number_ICs,
                                 sigma_noise)
  
  starting_pts_without_ICs <- rep(starting_pts,number_ICs)
  starting_pts_with_ICs <- starting_pts_without_ICs+noise_vector
  
  
  
  
  #this ensures any noisy initial values below 0 are changed to 0.01
  #similarly above 1 changed to 0.99
  starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs<0,0.01)
  starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs>1,0.99)
  
  
  
  
  #now we find the optim values for each alpha for each IC
  #we save the parameters and the optimised value
  optim_par <- rep(0,length(alphas)*number_ICs*batch_size*2)
  optim_value <- rep(0,length(alphas)*number_ICs)
  z <- 1
  for (i in 1:length(alphas)){
    for (j in 1:number_ICs){
      optim <- optim(method = "L-BFGS-B",
                     
                     par=starting_pts_with_ICs[(2*(j-1)*batch_size+1):(2*j*batch_size)], 
                     
                     fn=new_optim_fixed_runs,
                     
                     oldpts=batch0,
                     grid_length=grid_length,
                     theta_val2=theta_emulator,
                     sigma_val2=sigma_emulator,
                     E_f_val2=E_emulator,
                     nugget_val2=nugget_emulator,
                     true_runs_vec=batch0_evaluated,
                     alpha=alphas[i],
                     
                     runtime_emulator_theta=run_time_emulator_theta,
                     runtime_emulator_sigma=run_time_emulator_sigma,
                     runtime_emulator_E=run_time_emulator_E,
                     runtime_emulator_nugget=run_time_emulator_nugget,
                     
                     lower=0,
                     upper=1)
      optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)] <- optim$par
      optim_value[z] <- optim$value
      z <- z+1
    }
  }
  
  
  
  
  #this produces var graph showing locations of optimised points for every IC for every alpha
  #warning - this can produce hundreds of graphs!
  if (want_var_graph==TRUE){
    for (z in 1:(length(alphas)*number_ICs)){
      var_graph(xD_vec_original=batch0,
                xD_vec_new=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                grid_length=grid_length, 
                theta_val=theta_emulator, 
                sigma_val=sigma_emulator, 
                E_f_val=E_emulator,
                nugget_val=nugget_emulator)
    }
  }
  
  
  
  
  #we calculate v90 and Tc for each set of optimised locations
  if (want_grid==TRUE | want_graph==TRUE){
    #calculate V90 and Tc for each set of locations
    oldpts_matrix <- matrix(batch0,ncol=2,byrow=TRUE) 
    v90 <- rep(0,length(alphas)*number_ICs)
    Tc <- rep(0,length(alphas)*number_ICs)
    for (z in 1:(length(alphas)*number_ICs)){ #look at this section again... particularly the v90
      newpts_matrix <- matrix(optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],ncol=2,byrow=TRUE)
      xD <- rbind(oldpts_matrix,newpts_matrix)
      D <- f(xD)
      x_grid <- seq(-0.001,1.001,len=grid_length)
      xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
      em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_emulator,sigma=sigma_emulator,E_f=E_emulator,nugget=nugget_emulator)
      v90[z] <- quantile(em_out[,"VarD_f(x)"],percentile)
      Tc[z] <- emulate_log_run_time_Tc(new_batch_vector=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                                       xD_matrix=oldpts_matrix,
                                       D_vector=batch0_evaluated,
                                       theta_val=run_time_emulator_theta,
                                       sigma_val=run_time_emulator_sigma,
                                       E_f_val=run_time_emulator_E,
                                       nugget_val=run_time_emulator_nugget) #take care this should use the log emulator characterisitics!!! May have to change in future
    }
    
    
  }
  
  
  
  
  #for each alpha we find best locations (over the ICs) based on which gives the lowest optim value
  minimums_indices <- rep(0,length(alphas))
  for (k in 1:length(alphas)){
    minimums_indices[k] <- (k-1)*number_ICs+which.min(optim_value[((k-1)*number_ICs+1):(k*number_ICs)])
  }
  #pars_for_grid <- list_optim_par[[minimums_indices]] #this of length length(alphas)
  v90_for_grid <- v90[minimums_indices]
  Tc_for_grid <- Tc[minimums_indices]
  
  
  
  
  #a grid 
  #for each row e.g. 1/2 evaluate the relative criterion log10U = 1/2*log10V_{90}+1/2*log10T_{c} for the optim design using each alpha
  #if optim is performing/behaving correctly the diagonals should be the minimum of each row 
  if (want_grid==TRUE){
    #minimum for each different starting conditions
    matrix <- matrix(0,ncol=length(alphas),nrow=length(alphas))
    for (i in 1:length(alphas)){
      for (j in 1:length(alphas)){
        matrix[i,j] <- 1/2*log10(v90_for_grid[j])+alphas[i]*log10(Tc_for_grid[j])
      }
    }
    colnames(matrix) <- alphas
    rownames(matrix) <- alphas
    print(as.table(matrix))
  }
  
  
  
  
  #plot Tc on x axis against V90 on y axis for each alpha
  if (want_graph==TRUE){
    #graph
    par(mfrow=c(1,1))
    plot(Tc,v90,
         xlab="Tc",
         ylab="V90")
    #label the plot with fewer points
    plot(Tc_for_grid,v90_for_grid,
         xlab="Tc",
         ylab="V90")
    text(Tc_for_grid,v90_for_grid,label=round(alphas,2),cex= 0.7,pos=3)
  }
  
  
  
  
  #plots the variance graphs (showing where optimised locations are) only for the best locations for each alpha
  #alpha given in bottom left corner
  if (want_best_var_graphs==TRUE){
    for (z in 1:length(alphas)){
      var_graph(xD_vec_original=batch0,
                xD_vec_new=optim_par[((batch_size*2)*(minimums_indices[z]-1)+1):((batch_size*2)*minimums_indices[z])],
                grid_length=grid_length, 
                theta_val=theta_emulator, 
                sigma_val=sigma_emulator, 
                E_f_val=E_emulator,
                nugget_val=nugget_emulator)
      title(paste("Alpha =",round(alphas[z],2)),adj=0,line=-26)
    }
  }
  
  
  
  
  #return some things we may want to look at - especially the table (easy to save this way)
  return(list("Tc"=Tc,"v90"=v90,"Tc_for_grid"=Tc_for_grid,"v90_for_grid"=v90_for_grid,"minimums_indices"=minimums_indices,"optim_par"=optim_par,"optim_value"=optim_value,"table"=as.table(matrix)))
}



## end of run before alpha. ------------------------------------------------



## 8 pts -------------------------------------------------------------------

#png(file="big_function_8pts_%d.png",width=1400,height=1200,res=200)

big_function_8pts <- var_grid_graph(want_var_graph=FALSE,
                                    want_grid=TRUE,
                                    want_graph=TRUE,
                                    want_best_var_graphs=TRUE,
                                    
                                    batch0=c(t(xD)),
                                    batch0_evaluated=D,
                                    
                                    grid_length=60,
                                    
                                    batch_size=8, #4, 6, 8 or 10
                                    
                                    theta_emulator=1,
                                    sigma_emulator=sd(D),
                                    E_emulator=mean(D),
                                    nugget_emulator=0.145064^2,
                                    
                                    run_time_emulator_theta=1,
                                    run_time_emulator_sigma=sd(D),
                                    run_time_emulator_E=mean(D),
                                    run_time_emulator_nugget=0.145064^2,
                                    
                                    alphas=c(0,0.1,1/3,1/2,2/3,0.9,1,1.5,2,3,4), #vector of alphas we want
                                    percentile=0.9, #percentile of variance used in optim 
                                    number_ICs=10, #how many different starting points we want
                                    sigma_noise=0.1)


#dev.off()



write.csv(big_function_8pts$"table",file="/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/My work/RStudio/week12/big_function_8pts_table.csv")






# Figure 42 ---------------------------------------------------------------

## run before ceiling... ---------------------------------------------------
rm(list = ls())


library(viridisLite)
library(pdist)

nugget_argument_efficient_nugget_all_graphs <- function(xD_matrix, 
                                                        D_vector,
                                                        grid_length=50, 
                                                        theta_val=0.5, 
                                                        sigma_val=0.5, 
                                                        E_f_val=0,
                                                        
                                                        #new part:
                                                        nugget_val=10^(-6)){
  
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
  #new part (adding nugget as argument):
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)  
  
  
  #extracting the emulator expectation and variance
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  
  #graphing
  emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")
}





nugget_argument_efficient_nugget_simple_BL_emulator <- function(xP,             # the set of emulator prediction points
                                                                xD,             # the run input locations xD
                                                                D,              # the run outputs D = (f(x^1),...,f(x^n))
                                                                theta = 1,      # the correlation lengths (can be a vector)
                                                                sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                                                E_f = 0,        # prior expectation of f: E(f(x)) = 0 
                                                                using_pdist = 1,# if you have installed pdist package
                                                                
                                                                #new part:
                                                                nugget = 10^(-6) #must input nugget 0.1^2 
                                                                #squared because working with sigma squared
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
  
  
  
  
  
  
  #deterministic: manage numerical error
  #diagonal <- diag(10^-6,n,n) #1 millionth of what we are doing but well above numerical error (nice safe amount)
  #stochastic: manage stochasticity
  #diagonal <- diag(0.1^2,n,n) #each green pt quite a bit more uncertain
  
  #new part:
  diagonal <- diag(nugget,n,n)
  Var_D <- Var_D + sigma^2*diagonal
  
  #note that nugget 0.1^2 or 0.2^2 closer to 10% of mean of D (which is roughly what we want)
  #could put in directly but will have to change by hand every time change emulator parameters
  
  
  
  
  
  
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



dataframe_sim_run_times <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/simulation_overview.txt")
matrix_sim_run_times <- as.matrix(dataframe_sim_run_times)
colnames(matrix_sim_run_times) <- c("Emulator Coordinate x", 
                                    "Emulator Coordinate y", 
                                    "Runtime [CPU-Hours]", 
                                    "Batch Number", 
                                    "fmin",	
                                    "fmax",	
                                    "rho0 [cm^-3]",
                                    "T_AGN [K]")
x <- matrix_sim_run_times[,"Emulator Coordinate x"]
y <- matrix_sim_run_times[,"Emulator Coordinate y"]
plot(x,y)

xD <- matrix(data=c(x,y),ncol=2,byrow = FALSE)
D <- log10(matrix_sim_run_times[,"Runtime [CPU-Hours]"])

emulate_log_run_time_Tc <- function(new_batch_vector,
                                    xD_matrix,
                                    D_vector,
                                    
                                    theta_val,
                                    sigma_val,
                                    E_f_val,
                                    
                                    nugget_val=10^(-6))
  
{
  xD <- xD_matrix
  D <- D_vector
  new_batch_matrix <- matrix(data=new_batch_vector,ncol=2,byrow = TRUE)
  xP <- new_batch_matrix
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)
  E_D_fx_mat <- em_out[,"ExpD_f(x)"] 
  normal_times <- 10^(E_D_fx_mat)
  return(sum(normal_times))
}


var_graph <- function(xD_vec_original,
                      xD_vec_new,
                      grid_length=50, 
                      theta_val=0.5, 
                      sigma_val=0.5, 
                      E_f_val=0,
                      
                      nugget_val=10^(-6))
{
  par(mar=c(4,4,3,3))
  var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
  
  xD_vec <- c(xD_vec_original,xD_vec_new)
  xD <- matrix(xD_vec,ncol=2,byrow=TRUE)
  D <- f(xD)
  
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)   
  
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]",
                 #new part:
                 xD_col=rep(c("purple","pink"),c(length(xD_vec_original)/2,length(xD_vec_new)/2)) #/2 because x and y coords in vector so number of points is length/2
  )
  
}

#create random function (won't affect variance but function needs it)
f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))





library(MASS)

noise_function <- function(batch_size,
                           number_ICs,
                           sigma_noise){
  var_matrix <- matrix(c(sigma_noise^2,0,
                         0,sigma_noise^2),
                       nrow=2,
                       ncol=2,
                       byrow=TRUE
  )
  noise2D_matrix <- mvrnorm(n=batch_size*(number_ICs-1),
                            mu=rep(0,2),
                            Sigma=var_matrix)
  noise2D_vector <- c(rep(0,batch_size*2),c(t(noise2D_matrix)))
  return(noise2D_vector)
}






sigmoidal_optim_function <- function(newpts, #will optimise first argument
                                     oldpts,
                                     grid_length=50,
                                     theta_val2=0.5,
                                     sigma_val2=0.5,
                                     E_f_val2=0,
                                     nugget_val2=10^(-6),
                                     true_runs_vec,
                                     
                                     ceiling,
                                     freeze,
                                     
                                     runtime_emulator_theta,
                                     runtime_emulator_sigma,
                                     runtime_emulator_E,
                                     runtime_emulator_nugget
){
  oldpts_matrix <- matrix(oldpts,ncol=2,byrow=TRUE) 
  newpts_matrix <- matrix(newpts,ncol=2,byrow=TRUE)
  xD <- rbind(oldpts_matrix,newpts_matrix)
  D <- f(xD)
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val2,sigma=sigma_val2,E_f=E_f_val2,nugget=nugget_val2)
  v90 <- quantile(em_out[,"VarD_f(x)"],0.9)
  Tc <- emulate_log_run_time_Tc(new_batch_vector=newpts,
                                xD_matrix=oldpts_matrix,
                                D_vector=true_runs_vec,
                                theta_val=runtime_emulator_theta,
                                sigma_val=runtime_emulator_sigma,
                                E_f_val=runtime_emulator_E,
                                nugget_val=runtime_emulator_nugget)
  
  if (freeze==FALSE){
    a <- 5 #steepness
    b <- 3 #height
    c <- ceiling #log10(ceiling) x value halfway up the sigmoidal
    d <- 0.5 #steepness of linear term added on 
  }
  
  if (freeze==TRUE){
    a <- 10 #steepness - THIS IS THE ONLY DIFFERENCE
    b <- 3 #height
    c <- ceiling #log10(ceiling) x value halfway up the sigmoidal
    d <- 0.5 #steepness of linear term added on
  }
  
  x <- Tc
  
  if (Tc<=ceiling){
    sigmoidal <- b/(1+exp(-a*(x-c)))
  }
  if (Tc>ceiling){
    sigmoidal <- b/(1+exp(-a*(x-c))) + d*x - d*c
  }
  
  final <- (1/2)*log10(v90) + sigmoidal
  return(final)
}

# #checking sigmoidal part works ok
# sig <- function(log10Tc,
#                 log10ceiling
#                 ){
#   a <- 10 #5 or 10
#   b <- 3
#   c <- log10ceiling
#   d <- 0.5
#   x <- log10Tc
#   if (log10Tc<=log10ceiling){
#     sigmoidal <- b/(1+exp(-a*(x-c)))
#   }
#   if (log10Tc>log10ceiling){
#     sigmoidal <- b/(1+exp(-a*(x-c))) + d*x - d*c
#   }
#   return(sigmoidal)
# }
# x <- seq(0,10,by=0.01)
# which(x==4)
# plot(x,sapply(x,FUN=sig,log10ceiling=7),ylim=c(0,6))
# y <- sapply(x,FUN=sig,log10ceiling=7)
# y[601]

#Function using sigmoidal criteria
everything_function_sigmoidal <- function(want_all_var_graphs,
                                          want_only_best_var_graph,
                                          want_freeze, #initialise design at melt design found now with harsher sigmoidal
                                          
                                          batch0,
                                          batch0_evaluated, #evaluated run times
                                          
                                          grid_length,
                                          
                                          batch_size, #4, 6, 8 or 10 now 5,7 also allowed
                                          
                                          theta_emulator,
                                          sigma_emulator,
                                          E_emulator,
                                          nugget_emulator, #input as squared
                                          
                                          run_time_emulator_theta,
                                          run_time_emulator_sigma,
                                          run_time_emulator_E,
                                          run_time_emulator_nugget,
                                          
                                          percentile, #percentile of variance used in optim 
                                          number_ICs, #how many different starting points we want
                                          sigma_noise, #N(0,sigma_noise^2) for random noise
                                          
                                          ceiling_value #NOT LOGGED max total run time cannot exceed
){
  
  #create random function (won't affect variance but function needs it)
  f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))
  
  
  
  #determine initial starting points 
  #dependent on batch_size (number of new points) we want
  if (batch_size==4){
    starting_pts <- c(0.9,0.85,
                      0.6,0.95, 
                      0.4,0.05,
                      0.1,0.35)
  }
  if (batch_size==6){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05)
  }
  if (batch_size==8){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6)
  }
  if (batch_size==10){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6,
                      
                      0.2,0.8,
                      0.8,0.1)
  }
  if (batch_size==5){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.4,0.05,
                      0.1,0.35)
  }
  if (batch_size==7){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.7,0.6)
    
  }
  if (batch_size==9){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6,
                      
                      0.2,0.8)
  }
  
  
  
  
  #this creates "number_ICs" lots of initial starting points
  #1st is with no noise
  #2nd onwards are with 2D multivariate normal noise mean (0,0) var=(sigma_noise^2 , 0
  #                                                                  0 , sigma_noise^2)
  noise_vector <- noise_function(batch_size,
                                 number_ICs,
                                 sigma_noise)
  
  starting_pts_without_ICs <- rep(starting_pts,number_ICs)
  starting_pts_with_ICs <- starting_pts_without_ICs+noise_vector
  
  
  
  
  #this ensures any noisy initial values below 0 are changed to 0.01
  #similarly above 1 changed to 0.99
  starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs<0,0.01)
  starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs>1,0.99)
  
  
  
  
  #now we find the optim value over all ICs
  #we save the parameters and the optimised value
  optim_par <- rep(0,number_ICs*batch_size*2)
  optim_value <- rep(0,number_ICs)
  for (j in 1:number_ICs){
    optim <- optim(method = "L-BFGS-B",
                   
                   par=starting_pts_with_ICs[(2*(j-1)*batch_size+1):(2*j*batch_size)], 
                   
                   fn=sigmoidal_optim_function,
                   
                   oldpts=batch0,
                   grid_length=grid_length,
                   theta_val2=theta_emulator,
                   sigma_val2=sigma_emulator,
                   E_f_val2=E_emulator,
                   nugget_val2=nugget_emulator,
                   true_runs_vec=batch0_evaluated,
                   
                   runtime_emulator_theta=run_time_emulator_theta,
                   runtime_emulator_sigma=run_time_emulator_sigma,
                   runtime_emulator_E=run_time_emulator_E,
                   runtime_emulator_nugget=run_time_emulator_nugget,
                   
                   ceiling=ceiling_value,
                   freeze=FALSE, #this step is just about doing the melt (gentle) sigmoidal
                   
                   lower=0,
                   upper=1)
    optim_par[((batch_size*2)*(j-1)+1):((batch_size*2)*j)] <- optim$par
    optim_value[j] <- optim$value
  }
  
  
  
  
  
  
  #this produces var graph showing locations of optimised points for every IC for every alpha
  #warning - this can produce hundreds of graphs!
  if (want_all_var_graphs==TRUE){
    for (z in 1:(number_ICs)){
      var_graph(xD_vec_original=batch0,
                xD_vec_new=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                grid_length=grid_length, 
                theta_val=theta_emulator, 
                sigma_val=sigma_emulator, 
                E_f_val=E_emulator,
                nugget_val=nugget_emulator)
    }
  }
  
  
  
  
  #we calculate v90 and Tc for each set of optimised locations
  #calculate V90 and Tc for each set of locations
  oldpts_matrix <- matrix(batch0,ncol=2,byrow=TRUE) 
  v90 <- rep(0,number_ICs)
  Tc <- rep(0,number_ICs)
  for (z in 1:(number_ICs)){ 
    newpts_matrix <- matrix(optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],ncol=2,byrow=TRUE)
    xD <- rbind(oldpts_matrix,newpts_matrix)
    D <- f(xD)
    x_grid <- seq(-0.001,1.001,len=grid_length)
    xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
    em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_emulator,sigma=sigma_emulator,E_f=E_emulator,nugget=nugget_emulator)
    v90[z] <- quantile(em_out[,"VarD_f(x)"],percentile)
    Tc[z] <- emulate_log_run_time_Tc(new_batch_vector=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                                     xD_matrix=oldpts_matrix,
                                     D_vector=batch0_evaluated,
                                     theta_val=run_time_emulator_theta,
                                     sigma_val=run_time_emulator_sigma,
                                     E_f_val=run_time_emulator_E,
                                     nugget_val=run_time_emulator_nugget) #take care this should use the log emulator characterisitics!!! May have to change in future
  }
  
  
  
  
  
  #best design (over the ICs) based on which gives the lowest optim value
  minimum_index <- which.min(optim_value)
  
  best_design_v90 <- v90[minimum_index]
  best_design_Tc <- Tc[minimum_index]
  best_design_locations <- matrix(optim_par[((batch_size*2)*(minimum_index-1)+1):((batch_size*2)*minimum_index)],ncol=2,byrow=TRUE)
  best_design_value <- optim_value[minimum_index]
  
  
  
  #plots the variance graph (showing where optimised locations are) only for the best design
  if (want_only_best_var_graph==TRUE){
    var_graph(xD_vec_original=batch0,
              xD_vec_new=c(t(best_design_locations)),
              grid_length=grid_length, 
              theta_val=theta_emulator, 
              sigma_val=sigma_emulator, 
              E_f_val=E_emulator,
              nugget_val=nugget_emulator)
    title("Best Design",adj=0,line=-26)
    
  }
  
  
  
  
  freeze_v90 <- NA
  freeze_Tc <- NA
  freeze_locations <- rep(NA,batch_size)
  freeze_value <- NA
  
  if (want_freeze==TRUE){
    freeze_optim <- optim(method = "L-BFGS-B",
                          
                          par=c(t(best_design_locations)), 
                          
                          fn=sigmoidal_optim_function,
                          
                          oldpts=batch0,
                          grid_length=grid_length,
                          theta_val2=theta_emulator,
                          sigma_val2=sigma_emulator,
                          E_f_val2=E_emulator,
                          nugget_val2=nugget_emulator,
                          true_runs_vec=batch0_evaluated,
                          
                          runtime_emulator_theta=run_time_emulator_theta,
                          runtime_emulator_sigma=run_time_emulator_sigma,
                          runtime_emulator_E=run_time_emulator_E,
                          runtime_emulator_nugget=run_time_emulator_nugget,
                          
                          ceiling=ceiling_value,
                          freeze=TRUE, #now we use a stricter sigmoidal!
                          
                          lower=0,
                          upper=1)
    freeze_locations <- freeze_optim$par
    freeze_value <- freeze_optim$value
    
    var_graph(xD_vec_original=batch0,
              xD_vec_new=freeze_locations,
              grid_length=grid_length, 
              theta_val=theta_emulator, 
              sigma_val=sigma_emulator, 
              E_f_val=E_emulator,
              nugget_val=nugget_emulator)
    title("Freeze Design",adj=0,line=-26)
    
    oldpts_matrix <- matrix(batch0,ncol=2,byrow=TRUE)
    freeze_locations_matrix <- matrix(freeze_locations,ncol=2,byrow=TRUE)
    xD <- rbind(oldpts_matrix,freeze_locations_matrix)
    D <- f(xD)
    x_grid <- seq(-0.001,1.001,len=grid_length)
    xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
    em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_emulator,sigma=sigma_emulator,E_f=E_emulator,nugget=nugget_emulator)
    freeze_v90 <- quantile(em_out[,"VarD_f(x)"],percentile)
    freeze_Tc <- emulate_log_run_time_Tc(new_batch_vector=freeze_locations,
                                         xD_matrix=oldpts_matrix,
                                         D_vector=batch0_evaluated,
                                         theta_val=run_time_emulator_theta,
                                         sigma_val=run_time_emulator_sigma,
                                         E_f_val=run_time_emulator_E,
                                         nugget_val=run_time_emulator_nugget)
  }
  
  
  
  #return some things we may want to look at
  return(list("All v90"=v90,
              "All Tc"=Tc,
              "All locations"=optim_par,
              "All values"=optim_value,
              
              "Best Design v90"=best_design_v90,
              "Best Design Tc"=best_design_Tc,
              "Best Design locations"=best_design_locations,
              "Best Design value"=best_design_value,
              
              "Freeze v90"=freeze_v90,
              "Freeze Tc"=freeze_Tc,
              "Freeze location"=freeze_locations,
              "Freeze value"=freeze_value,
              
              "Best design index"=minimum_index
  ))
}


## end of run before ceiling. ----------------------------------------------



## run 6 points (theta=0.5) 100 ICs hist v90 hist Tc -----------------------
#png(file="q_ceiling_100ics_%d.png",width=1400,height=1200,res=200)
q_ceiling_100ics <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                  want_only_best_var_graph=TRUE,
                                                  want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                  
                                                  batch0=c(t(xD)),
                                                  batch0_evaluated=D, #evaluated run times
                                                  
                                                  grid_length=60,
                                                  
                                                  batch_size=6, #4, 6, 8 or 10
                                                  
                                                  theta_emulator=1/2,
                                                  sigma_emulator=sd(D),
                                                  E_emulator=mean(D),
                                                  nugget_emulator=0.001^2, #input as squared
                                                  
                                                  run_time_emulator_theta=1,
                                                  run_time_emulator_sigma=sd(D),
                                                  run_time_emulator_E=mean(D),
                                                  run_time_emulator_nugget=0.145064^2,
                                                  
                                                  percentile=0.9, #percentile of variance used in optim 
                                                  number_ICs=100, #how many different starting points we want
                                                  sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                  
                                                  ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()

#png(file="q_ceiling_histograms_100ics_%d.png",width=1400,height=1200,res=200)
hist(q_ceiling_100ics$"All v90",
     xlab="V90",
     ylab="Counts",
     main = "Histogram of V90")
hist(q_ceiling_100ics$"All Tc",
     xlab="Tc",
     ylab="Counts",
     main="Histogram of Tc")
#dev.off()
q_ceiling_100ics$"Best Design v90" #0.0004790062
q_ceiling_100ics$"Best Design Tc" #5924.791






## smaller bw: run 6 points (theta=0.5) 100 ICs hist v90 hist Tc -----------------------
#png(file="smallerbw_q_ceiling_100ics_%d.png",width=1400,height=1200,res=200)
q_ceiling_100ics <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                                  want_only_best_var_graph=TRUE,
                                                  want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                                  
                                                  batch0=c(t(xD)),
                                                  batch0_evaluated=D, #evaluated run times
                                                  
                                                  grid_length=60,
                                                  
                                                  batch_size=6, #4, 6, 8 or 10
                                                  
                                                  theta_emulator=1/2,
                                                  sigma_emulator=sd(D),
                                                  E_emulator=mean(D),
                                                  nugget_emulator=0.001^2, #input as squared
                                                  
                                                  run_time_emulator_theta=1,
                                                  run_time_emulator_sigma=sd(D),
                                                  run_time_emulator_E=mean(D),
                                                  run_time_emulator_nugget=0.145064^2,
                                                  
                                                  percentile=0.9, #percentile of variance used in optim 
                                                  number_ICs=100, #how many different starting points we want
                                                  sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                                  
                                                  ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)
#dev.off()

a <- q_ceiling_100ics$"All v90"
b <- q_ceiling_100ics$"All Tc"
c <- q_ceiling_100ics$"All locations"
d <- q_ceiling_100ics$"All values"
e <- q_ceiling_100ics$"Best Design v90"
f <- q_ceiling_100ics$"Best Design Tc"
g <- q_ceiling_100ics$"Best Design locations"
h <- q_ceiling_100ics$"Best Design value"
i <- q_ceiling_100ics$"Freeze v90"
j <- q_ceiling_100ics$"Freeze Tc"
k <- q_ceiling_100ics$"Freeze location"
l <- q_ceiling_100ics$"Freeze value"
m <- q_ceiling_100ics$"Best design index"

save(a, 
     b, 
     c, 
     d, 
     e, 
     f, 
     g, 
     h, 
     i, 
     j, 
     k, 
     l, 
     m, 
     file = "q_ceiling_100ics.RData")

load("q_ceiling_100ics.RData")

#png(file="other_option_smaller_bw_q_ceiling_histograms_100ics_%d.png",width=1400,height=1200,res=200)
hist(a,
     xlab="V90",
     ylab="Counts",
     main = "Histogram of V90",
     breaks=seq(0,0.006,by=0.0004),
     xlim=c(0,0.006))
#dev.off()



#png(file="smaller_bw_q_ceiling_histograms_100ics_%d.png",width=1400,height=1200,res=200)
hist(q_ceiling_100ics$"All v90",
     xlab="V90",
     ylab="Counts",
     main = "Histogram of V90",
     breaks=seq(0,0.006,by=0.0005))
hist(q_ceiling_100ics$"All Tc",
     xlab="Tc",
     ylab="Counts",
     main="Histogram of Tc",
     breaks=seq(3000,6000,by=250))
plot(q_ceiling_100ics$"All Tc",
     q_ceiling_100ics$"All v90",
     xlab="Tc",
     ylab="V90")
#dev.off()
q_ceiling_100ics$"Best Design v90" #0.000466063
q_ceiling_100ics$"Best Design Tc" #5909.531




#png(file="larger_smaller_bw_q_ceiling_histograms_100ics_%d.png",width=1400,height=1200,res=200)
hist(q_ceiling_100ics$"All v90",
     xlab="V90",
     ylab="Counts",
     main = "Histogram of V90",
     breaks=seq(0,0.006,by=0.00075))
#dev.off()








# Figure 43 ---------------------------------------------------------------

## 43a ---------------------------------------------------------------------

# See the 6-point design with 6000 ceiling
# Figure 23c


## 43b ---------------------------------------------------------------------


### run all functions needed --------------------------------------------------------------------



rm(list = ls())


library(viridisLite)
library(pdist)

nugget_argument_efficient_nugget_all_graphs <- function(xD_matrix, 
                                                        D_vector,
                                                        grid_length=50, 
                                                        theta_val=0.5, 
                                                        sigma_val=0.5, 
                                                        E_f_val=0,
                                                        
                                                        #new part:
                                                        nugget_val=10^(-6)){
  
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
  #new part (adding nugget as argument):
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)  
  
  
  #extracting the emulator expectation and variance
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  
  #graphing
  emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=exp_cols,main="Emulator Adjusted Expectation E_D[f(x)]")
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]")
}





nugget_argument_efficient_nugget_simple_BL_emulator <- function(xP,             # the set of emulator prediction points
                                                                xD,             # the run input locations xD
                                                                D,              # the run outputs D = (f(x^1),...,f(x^n))
                                                                theta = 1,      # the correlation lengths (can be a vector)
                                                                sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                                                E_f = 0,        # prior expectation of f: E(f(x)) = 0 
                                                                using_pdist = 1,# if you have installed pdist package
                                                                
                                                                #new part:
                                                                nugget = 10^(-6) #must input nugget 0.1^2 
                                                                #squared because working with sigma squared
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
  
  
  
  
  
  
  #deterministic: manage numerical error
  #diagonal <- diag(10^-6,n,n) #1 millionth of what we are doing but well above numerical error (nice safe amount)
  #stochastic: manage stochasticity
  #diagonal <- diag(0.1^2,n,n) #each green pt quite a bit more uncertain
  
  #new part:
  diagonal <- diag(nugget,n,n)
  Var_D <- Var_D + sigma^2*diagonal
  
  #note that nugget 0.1^2 or 0.2^2 closer to 10% of mean of D (which is roughly what we want)
  #could put in directly but will have to change by hand every time change emulator parameters
  
  
  
  
  
  
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



dataframe_sim_run_times <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/simulation_overview.txt")
matrix_sim_run_times <- as.matrix(dataframe_sim_run_times)
colnames(matrix_sim_run_times) <- c("Emulator Coordinate x", 
                                    "Emulator Coordinate y", 
                                    "Runtime [CPU-Hours]", 
                                    "Batch Number", 
                                    "fmin",	
                                    "fmax",	
                                    "rho0 [cm^-3]",
                                    "T_AGN [K]")
x <- matrix_sim_run_times[,"Emulator Coordinate x"]
y <- matrix_sim_run_times[,"Emulator Coordinate y"]
plot(x,y)

xD <- matrix(data=c(x,y),ncol=2,byrow = FALSE)
D <- log10(matrix_sim_run_times[,"Runtime [CPU-Hours]"])

emulate_log_run_time_Tc <- function(new_batch_vector,
                                    xD_matrix,
                                    D_vector,
                                    
                                    theta_val,
                                    sigma_val,
                                    E_f_val,
                                    
                                    nugget_val=10^(-6))
  
{
  xD <- xD_matrix
  D <- D_vector
  new_batch_matrix <- matrix(data=new_batch_vector,ncol=2,byrow = TRUE)
  xP <- new_batch_matrix
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)
  E_D_fx_mat <- em_out[,"ExpD_f(x)"] 
  normal_times <- 10^(E_D_fx_mat)
  return(sum(normal_times))
}


var_graph <- function(xD_vec_original,
                      xD_vec_new,
                      grid_length=50, 
                      theta_val=0.5, 
                      sigma_val=0.5, 
                      E_f_val=0,
                      
                      nugget_val=10^(-6))
{
  par(mar=c(4,4,3,3))
  var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
  
  xD_vec <- c(xD_vec_original,xD_vec_new)
  xD <- matrix(xD_vec,ncol=2,byrow=TRUE)
  D <- f(xD)
  
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)   
  
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]",
                 #new part:
                 xD_col=rep(c("purple","pink"),c(length(xD_vec_original)/2,length(xD_vec_new)/2)) #/2 because x and y coords in vector so number of points is length/2
  )
  
}

#create random function (won't affect variance but function needs it)
f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))





library(MASS)

noise_function <- function(batch_size,
                           number_ICs,
                           sigma_noise){
  var_matrix <- matrix(c(sigma_noise^2,0,
                         0,sigma_noise^2),
                       nrow=2,
                       ncol=2,
                       byrow=TRUE
  )
  noise2D_matrix <- mvrnorm(n=batch_size*(number_ICs-1),
                            mu=rep(0,2),
                            Sigma=var_matrix)
  noise2D_vector <- c(rep(0,batch_size*2),c(t(noise2D_matrix)))
  return(noise2D_vector)
}






simulation_Tc <- function(sim_new_pts_vec,
                          sim_old_pts_vec,
                          sim_runtime_D,
                          sim_runtime_theta,
                          sim_runtime_sigma,
                          sim_runtime_E,
                          sim_runtime_nugget,
                          sim_batchsize){
  
  xD <- matrix(sim_old_pts_vec,ncol=2,byrow = TRUE)
  D <- sim_runtime_D
  xP <- matrix(sim_new_pts_vec,ncol=2,byrow = TRUE)
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=sim_runtime_theta,sigma=sim_runtime_sigma,E_f=sim_runtime_E,nugget=sim_runtime_nugget)
  emulator_expectation <- em_out[,"ExpD_f(x)"]
  emulator_variance <- em_out[,"VarD_f(x)"]
  
  expected_total <- sum(10^emulator_expectation)
  
  simulations <- matrix(0,ncol=1000,nrow=sim_batchsize)
  for (i in 1:sim_batchsize){
    #MUST INPUT STANDARD DEVIATION!!!!!!!
    log_simulations <- rnorm(n=1000, mean=emulator_expectation[i], sd=sqrt(emulator_variance[i])) #MUST INPUT STANDARD DEVIATION!!!!!!!
    simulations[i,] <- 10^log_simulations
  }
  total_simulations <- colSums(simulations)
  mean_total <- mean(total_simulations)
  return(list("expected total"=expected_total, "simulated totals"=total_simulations, "mean total"=mean_total))
}

#I created this fn for hm but this will also be useful:
Tc_sigmoidal <- function(Tc_val,ceiling_val){
  a <- 5 #steepness
  b <- 3 #height
  c <- ceiling_val #where centred (x value halfway up the sigmoidal)
  d <- 0.5 #steepness of linear term added on
  x <- Tc_val
  if (Tc_val<=c){
    sigmoidal <- b/(1+exp(-a*(x-c)))
  }
  if (Tc_val>c){
    sigmoidal <- b/(1+exp(-a*(x-c))) + d*x - d*c
  }
  return(sigmoidal)
}

#now changing optim criteria to incorporate this function:
new_simulation_sigmoidal_optim_function <- function(newpts, #will optimise first argument
                                                    oldpts,
                                                    grid_length=50,
                                                    theta_val2=0.5,
                                                    sigma_val2=0.5,
                                                    E_f_val2=0,
                                                    nugget_val2=10^(-6),
                                                    true_runs_vec,
                                                    
                                                    ceiling,
                                                    freeze,
                                                    
                                                    runtime_emulator_theta,
                                                    runtime_emulator_sigma,
                                                    runtime_emulator_E,
                                                    runtime_emulator_nugget,
                                                    
                                                    #NEW!
                                                    percentile_val3,
                                                    
                                                    val3_batchsize
){
  oldpts_matrix <- matrix(oldpts,ncol=2,byrow=TRUE) 
  newpts_matrix <- matrix(newpts,ncol=2,byrow=TRUE)
  xD <- rbind(oldpts_matrix,newpts_matrix)
  D <- f(xD)
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val2,sigma=sigma_val2,E_f=E_f_val2,nugget=nugget_val2)
  v90 <- quantile(em_out[,"VarD_f(x)"],percentile_val3)
  
  # Tc <- emulate_log_run_time_Tc(new_batch_vector=newpts,
  #                               xD_matrix=oldpts_matrix,
  #                               D_vector=true_runs_vec,
  #                               theta_val=runtime_emulator_theta,
  #                               sigma_val=runtime_emulator_sigma,
  #                               E_f_val=runtime_emulator_E,
  #                               nugget_val=runtime_emulator_nugget)
  # 
  # if (freeze==FALSE){
  #   a <- 5 #steepness
  #   b <- 3 #height
  #   c <- ceiling #log10(ceiling) x value halfway up the sigmoidal
  #   d <- 0.5 #steepness of linear term added on 
  # }
  # 
  # if (freeze==TRUE){
  #   a <- 10 #steepness - THIS IS THE ONLY DIFFERENCE
  #   b <- 3 #height
  #   c <- ceiling #log10(ceiling) x value halfway up the sigmoidal
  #   d <- 0.5 #steepness of linear term added on
  # }
  # 
  # x <- Tc
  # 
  # if (Tc<=ceiling){
  #   sigmoidal <- b/(1+exp(-a*(x-c)))
  # }
  # if (Tc>ceiling){
  #   sigmoidal <- b/(1+exp(-a*(x-c))) + d*x - d*c
  # }
  
  #need to change the line section above
  
  Tc_function <- simulation_Tc(sim_new_pts_vec=newpts,
                               sim_old_pts_vec=oldpts,
                               sim_runtime_D=true_runs_vec,
                               sim_runtime_theta=runtime_emulator_theta,
                               sim_runtime_sigma=runtime_emulator_sigma,
                               sim_runtime_E=runtime_emulator_E,
                               sim_runtime_nugget=runtime_emulator_nugget,
                               sim_batchsize=val3_batchsize)
  Tc_simulations <- Tc_function$"simulated totals"
  sigmoidal_Tc_simulations <- sapply(Tc_simulations, Tc_sigmoidal, ceiling_val=ceiling)
  
  sigmoidal <- mean(sigmoidal_Tc_simulations)
  
  final <- (1/2)*log10(v90) + sigmoidal
  return(final)
}

new_simulation_everything_function_sigmoidal <- function(want_all_var_graphs,
                                                         want_only_best_var_graph,
                                                         want_freeze, #initialise design at melt design found now with harsher sigmoidal
                                                         
                                                         batch0,
                                                         batch0_evaluated, #evaluated run times
                                                         
                                                         grid_length,
                                                         
                                                         batch_size, #4, 6, 8 or 10 now 5,7 also allowed
                                                         
                                                         theta_emulator,
                                                         sigma_emulator,
                                                         E_emulator,
                                                         nugget_emulator, #input as squared
                                                         
                                                         run_time_emulator_theta,
                                                         run_time_emulator_sigma,
                                                         run_time_emulator_E,
                                                         run_time_emulator_nugget,
                                                         
                                                         percentile, #percentile of variance used in optim 
                                                         number_ICs, #how many different starting points we want
                                                         sigma_noise, #N(0,sigma_noise^2) for random noise
                                                         
                                                         ceiling_value #NOT LOGGED max total run time cannot exceed
){
  
  #create random function (won't affect variance but function needs it)
  f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))
  
  
  
  #determine initial starting points 
  #dependent on batch_size (number of new points) we want
  if (batch_size==4){
    starting_pts <- c(0.9,0.85,
                      0.6,0.95, 
                      0.4,0.05,
                      0.1,0.35)
  }
  if (batch_size==6){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05)
  }
  if (batch_size==8){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6)
  }
  if (batch_size==10){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6,
                      
                      0.2,0.8,
                      0.8,0.1)
  }
  if (batch_size==5){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.4,0.05,
                      0.1,0.35)
  }
  if (batch_size==7){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.7,0.6)
    
  }
  if (batch_size==9){
    starting_pts <- c(0.6,0.95,
                      0.9,0.85,
                      0.95,0.65,
                      
                      0.05,0.4,
                      0.1,0.15,
                      0.35,0.05,
                      
                      0.3,0.5,
                      0.7,0.6,
                      
                      0.2,0.8)
  }
  
  
  
  
  #this creates "number_ICs" lots of initial starting points
  #1st is with no noise
  #2nd onwards are with 2D multivariate normal noise mean (0,0) var=(sigma_noise^2 , 0
  #                                                                  0 , sigma_noise^2)
  noise_vector <- noise_function(batch_size,
                                 number_ICs,
                                 sigma_noise)
  
  starting_pts_without_ICs <- rep(starting_pts,number_ICs)
  starting_pts_with_ICs <- starting_pts_without_ICs+noise_vector
  
  
  
  
  #this ensures any noisy initial values below 0 are changed to 0.01
  #similarly above 1 changed to 0.99
  starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs<0,0.01)
  starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs>1,0.99)
  
  
  
  
  #now we find the optim value over all ICs
  #we save the parameters and the optimised value
  optim_par <- rep(0,number_ICs*batch_size*2)
  optim_value <- rep(0,number_ICs)
  for (j in 1:number_ICs){
    optim <- optim(method = "L-BFGS-B",
                   
                   par=starting_pts_with_ICs[(2*(j-1)*batch_size+1):(2*j*batch_size)], 
                   
                   fn=new_simulation_sigmoidal_optim_function,
                   
                   oldpts=batch0,
                   grid_length=grid_length,
                   theta_val2=theta_emulator,
                   sigma_val2=sigma_emulator,
                   E_f_val2=E_emulator,
                   nugget_val2=nugget_emulator,
                   true_runs_vec=batch0_evaluated,
                   
                   runtime_emulator_theta=run_time_emulator_theta,
                   runtime_emulator_sigma=run_time_emulator_sigma,
                   runtime_emulator_E=run_time_emulator_E,
                   runtime_emulator_nugget=run_time_emulator_nugget,
                   
                   ceiling=ceiling_value,
                   freeze=FALSE, #this step is just about doing the melt (gentle) sigmoidal
                   
                   percentile_val3 = percentile,
                   
                   val3_batchsize=batch_size,
                   
                   lower=0,
                   upper=1)
    optim_par[((batch_size*2)*(j-1)+1):((batch_size*2)*j)] <- optim$par
    optim_value[j] <- optim$value
  }
  
  
  
  
  
  
  #this produces var graph showing locations of optimised points for every IC for every alpha
  #warning - this can produce hundreds of graphs!
  if (want_all_var_graphs==TRUE){
    for (z in 1:(number_ICs)){
      var_graph(xD_vec_original=batch0,
                xD_vec_new=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                grid_length=grid_length, 
                theta_val=theta_emulator, 
                sigma_val=sigma_emulator, 
                E_f_val=E_emulator,
                nugget_val=nugget_emulator)
    }
  }
  
  
  
  
  #we calculate v90 and Tc for each set of optimised locations
  #calculate V90 and Tc for each set of locations
  oldpts_matrix <- matrix(batch0,ncol=2,byrow=TRUE) 
  v90 <- rep(0,number_ICs)
  Tc <- rep(0,number_ICs)
  for (z in 1:(number_ICs)){ 
    newpts_matrix <- matrix(optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],ncol=2,byrow=TRUE)
    xD <- rbind(oldpts_matrix,newpts_matrix)
    D <- f(xD)
    x_grid <- seq(-0.001,1.001,len=grid_length)
    xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
    em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_emulator,sigma=sigma_emulator,E_f=E_emulator,nugget=nugget_emulator)
    v90[z] <- quantile(em_out[,"VarD_f(x)"],percentile)
    Tc[z] <- emulate_log_run_time_Tc(new_batch_vector=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                                     xD_matrix=oldpts_matrix,
                                     D_vector=batch0_evaluated,
                                     theta_val=run_time_emulator_theta,
                                     sigma_val=run_time_emulator_sigma,
                                     E_f_val=run_time_emulator_E,
                                     nugget_val=run_time_emulator_nugget) #take care this should use the log emulator characterisitics!!! May have to change in future
  }
  
  
  
  
  
  #best design (over the ICs) based on which gives the lowest optim value
  minimum_index <- which.min(optim_value)
  
  best_design_v90 <- v90[minimum_index]
  best_design_Tc <- Tc[minimum_index]
  best_design_locations <- matrix(optim_par[((batch_size*2)*(minimum_index-1)+1):((batch_size*2)*minimum_index)],ncol=2,byrow=TRUE)
  best_design_value <- optim_value[minimum_index]
  
  
  
  #plots the variance graph (showing where optimised locations are) only for the best design
  if (want_only_best_var_graph==TRUE){
    var_graph(xD_vec_original=batch0,
              xD_vec_new=c(t(best_design_locations)),
              grid_length=grid_length, 
              theta_val=theta_emulator, 
              sigma_val=sigma_emulator, 
              E_f_val=E_emulator,
              nugget_val=nugget_emulator)
    title("Best Design",adj=0,line=-26)
    
  }
  
  
  
  
  freeze_v90 <- NA
  freeze_Tc <- NA
  freeze_locations <- rep(NA,batch_size)
  freeze_value <- NA
  
  if (want_freeze==TRUE){
    freeze_optim <- optim(method = "L-BFGS-B",
                          
                          par=c(t(best_design_locations)), 
                          
                          fn=sigmoidal_optim_function,
                          
                          oldpts=batch0,
                          grid_length=grid_length,
                          theta_val2=theta_emulator,
                          sigma_val2=sigma_emulator,
                          E_f_val2=E_emulator,
                          nugget_val2=nugget_emulator,
                          true_runs_vec=batch0_evaluated,
                          
                          runtime_emulator_theta=run_time_emulator_theta,
                          runtime_emulator_sigma=run_time_emulator_sigma,
                          runtime_emulator_E=run_time_emulator_E,
                          runtime_emulator_nugget=run_time_emulator_nugget,
                          
                          ceiling=ceiling_value,
                          freeze=TRUE, #now we use a stricter sigmoidal!
                          
                          lower=0,
                          upper=1)
    freeze_locations <- freeze_optim$par
    freeze_value <- freeze_optim$value
    
    var_graph(xD_vec_original=batch0,
              xD_vec_new=freeze_locations,
              grid_length=grid_length, 
              theta_val=theta_emulator, 
              sigma_val=sigma_emulator, 
              E_f_val=E_emulator,
              nugget_val=nugget_emulator)
    title("Freeze Design",adj=0,line=-26)
    
    oldpts_matrix <- matrix(batch0,ncol=2,byrow=TRUE)
    freeze_locations_matrix <- matrix(freeze_locations,ncol=2,byrow=TRUE)
    xD <- rbind(oldpts_matrix,freeze_locations_matrix)
    D <- f(xD)
    x_grid <- seq(-0.001,1.001,len=grid_length)
    xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
    em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_emulator,sigma=sigma_emulator,E_f=E_emulator,nugget=nugget_emulator)
    freeze_v90 <- quantile(em_out[,"VarD_f(x)"],percentile)
    freeze_Tc <- emulate_log_run_time_Tc(new_batch_vector=freeze_locations,
                                         xD_matrix=oldpts_matrix,
                                         D_vector=batch0_evaluated,
                                         theta_val=run_time_emulator_theta,
                                         sigma_val=run_time_emulator_sigma,
                                         E_f_val=run_time_emulator_E,
                                         nugget_val=run_time_emulator_nugget)
  }
  
  
  
  #return some things we may want to look at
  return(list("All v90"=v90,
              "All Tc"=Tc,
              "All locations"=optim_par,
              "All values"=optim_value,
              
              "Best Design v90"=best_design_v90,
              "Best Design Tc"=best_design_Tc,
              "Best Design locations"=best_design_locations,
              "Best Design value"=best_design_value,
              
              "Freeze v90"=freeze_v90,
              "Freeze Tc"=freeze_Tc,
              "Freeze location"=freeze_locations,
              "Freeze value"=freeze_value,
              
              "Best design index"=minimum_index
  ))
}





## design ------------------------------------------------------------------


#ceiling will kick in
#compare to 6000 6 pts found previously!

#set ceiling to 6500 (500 above like histogram suggests)
#now see if similar to before...

# png(file="simulations_pts6_6500_%d.png",width=1400,height=1200,res=200)
# simulations_pts6_6500 <- new_simulation_everything_function_sigmoidal(want_all_var_graphs=FALSE,
#                                                                       want_only_best_var_graph=TRUE,
#                                                                       want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
# 
#                                                                       batch0=c(t(xD)),
#                                                                       batch0_evaluated=D, #evaluated run times
# 
#                                                                       grid_length=60,
# 
#                                                                       batch_size=6, #4, 6, 8 or 10
# 
#                                                                       theta_emulator=1/2,
#                                                                       sigma_emulator=sd(D),
#                                                                       E_emulator=mean(D),
#                                                                       nugget_emulator=0.001^2, #input as squared
# 
#                                                                       run_time_emulator_theta=1,
#                                                                       run_time_emulator_sigma=sd(D),
#                                                                       run_time_emulator_E=mean(D),
#                                                                       run_time_emulator_nugget=0.145064^2,
# 
#                                                                       percentile=0.9, #percentile of variance used in optim
#                                                                       number_ICs=40, #how many different starting points we want
#                                                                       sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
# 
#                                                                       ceiling_value=6500 #NOT LOGGED max total run time cannot exceed
# )
# dev.off()
# saveRDS(simulations_pts6_6500, file="simulations_pts6_6500.RData")





# End of Appendices -------------------------------------------------------


