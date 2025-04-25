
# Section 7: Applications to Cosmological Simulations ---------------------

# Figures 21-24

# Functions ---------------------------------------------------------------

# Figure 21 ---------------------------------------------------------------

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

#png(file="batch0.png",width=1200,height=1315,res=200)
plot(x,y,
     xlim=c(0,1),
     ylim=c(0,1),
     pch = 19,
     main = "Initial Batch")
#dev.off()






# Figure 22 ---------------------------------------------------------------

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




xD <- matrix(data=c(x,y),ncol=2,byrow = FALSE)
D <- log10(matrix_sim_run_times[,"Runtime [CPU-Hours]"])
mean(D) #using this as our E_f_val
sd(D) #use this for our sigma

#which   ?*sd(log(D)) = sd(log data Shaun has sent) 
#file.choose()
stochasticity_data <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stochasticity_tests/Runtimes.txt")
log10_stochasticity_data <- log10(stochasticity_data$V1)
sd(log10_stochasticity_data) #0.02600337
sd(D) #0.1792545
sd(log10_stochasticity_data)/sd(D) #0.145064
#therefore use 0.145064^2 as our nugget term
#equivalently (sd(log10_stochasticity_data)^2)/(sd(D)^2)

#png(file="emulator_log_runtime_%d.png",width=1400,height=1200,res=200)
nugget_argument_efficient_nugget_all_graphs(xD_matrix = xD,
                                            D_vector = D,
                                            grid_length = 50,
                                            theta_val = 1, #pick a large theta BECAUSE WE EXPECT RUN TIME FUNCTION TO BE SMOOTH!
                                            sigma_val = sd(D), 
                                            E_f_val = mean(D),
                                            
                                            nugget_val = 0.145064^2) #must input it as squared as dealing with sigma squared
#dev.off()







# Run this before Figure 23 -----------------------------------------------

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











# Figure 23 ---------------------------------------------------------------

## 6000 ceiling 4 points ---------------------------------------------------

#png(file="pts4_6000_%d.png",width=1400,height=1200,res=200)
pts4_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                           want_only_best_var_graph=TRUE,
                                           want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                           
                                           batch0=c(t(xD)),
                                           batch0_evaluated=D, #evaluated run times
                                           
                                           grid_length=60,
                                           
                                           batch_size=4, #4, 6, 8 or 10
                                           
                                           theta_emulator=1/2,
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

# v90 <- pts4_6000$"All v90"
# Tc <- pts4_6000$"All Tc"
# values <- pts4_6000$"All values"
# hist(v90)
# hist(Tc)
# hist(values)
# plot(Tc,v90,
#      xlab="Tc",
#      ylab="V_{90}")

#dev.off()

#make a note of these!!!!!  
pts4_6000$"Best Design v90" #0.0008941561
pts4_6000$"Best Design Tc" #4236.612
pts4_6000$"Best Design value" #-1.524293



## 6000 ceiling 5 points ---------------------------------------------------

#png(file="pts5_6000_%d.png",width=1400,height=1200,res=200)
pts5_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                           want_only_best_var_graph=TRUE,
                                           want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                           
                                           batch0=c(t(xD)),
                                           batch0_evaluated=D, #evaluated run times
                                           
                                           grid_length=60,
                                           
                                           batch_size=5, #4, 6, 8 or 10
                                           
                                           theta_emulator=1/2,
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

# v90 <- pts5_6000$"All v90"
# Tc <- pts5_6000$"All Tc"
# values <- pts5_6000$"All values"
# hist(v90)
# hist(Tc)
# hist(values)
# plot(Tc,V90,
#      xlab="Tc",
#      ylab="V_{90}")

#dev.off()

#make a note of these!!!!!  
pts5_6000$"Best Design v90" #0.0006100015
pts5_6000$"Best Design Tc" #4784.707
pts5_6000$"Best Design value" #-1.607335






## 6000 ceiling 6 points ---------------------------------------------------

#png(file="pts6_6000_%d.png",width=1400,height=1200,res=200)
pts6_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
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
                                           number_ICs=40, #how many different starting points we want
                                           sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
                                           
                                           ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
)

# v90 <- pts6_6000$"All v90"
# Tc <- pts6_6000$"All Tc"
# values <- pts6_6000$"All values"
# hist(v90)
# hist(Tc)
# hist(values)
# plot(Tc,V90,
#      xlab="Tc",
#      ylab="V_{90}")

#dev.off()

#make a note of these!!!!!  
pts6_6000$"Best Design v90" #0.000475434
pts6_6000$"Best Design Tc" #5906.787
pts6_6000$"Best Design value" #-1.661455



## 6000 ceiling 7 points ---------------------------------------------------

#png(file="pts7_6000_%d.png",width=1400,height=1200,res=200)
pts7_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                           want_only_best_var_graph=TRUE,
                                           want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                           
                                           batch0=c(t(xD)),
                                           batch0_evaluated=D, #evaluated run times
                                           
                                           grid_length=60,
                                           
                                           batch_size=7, #4, 6, 8 or 10
                                           
                                           theta_emulator=1/2,
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

# v90 <- pts7_6000$"All v90"
# Tc <- pts7_6000$"All Tc"
# values <- pts7_6000$"All values"
# hist(v90)
# hist(Tc)
# hist(values)
# plot(Tc,V90,
#      xlab="Tc",
#      ylab="V_{90}")

#dev.off()

#make a note of these!!!!!  
pts7_6000$"Best Design v90" #0.0004916595
pts7_6000$"Best Design Tc" #5977.932
pts7_6000$"Best Design value" #-1.654168





## 6000 ceiling 8 points ---------------------------------------------------

#png(file="pts8_6000_%d.png",width=1400,height=1200,res=200)
pts8_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                           want_only_best_var_graph=TRUE,
                                           want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                           
                                           batch0=c(t(xD)),
                                           batch0_evaluated=D, #evaluated run times
                                           
                                           grid_length=60,
                                           
                                           batch_size=8, #4, 6, 8 or 10
                                           
                                           theta_emulator=1/2,
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

# v90 <- pts8_6000$"All v90"
# Tc <- pts8_6000$"All Tc"
# values <- pts8_6000$"All values"
# hist(v90)
# hist(Tc)
# hist(values)
# plot(Tc,V90,
#      xlab="Tc",
#      ylab="V_{90}")

#dev.off()

#make a note of these!!!!!  
pts8_6000$"Best Design v90" #0.0006746213
pts8_6000$"Best Design Tc" #5991.483
pts8_6000$"Best Design value" #-1.58547





## 6000 ceiling 9 points ---------------------------------------------------

#png(file="pts9_6000_%d.png",width=1400,height=1200,res=200)
pts9_6000 <- everything_function_sigmoidal(want_all_var_graphs=FALSE,
                                           want_only_best_var_graph=TRUE,
                                           want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
                                           
                                           batch0=c(t(xD)),
                                           batch0_evaluated=D, #evaluated run times
                                           
                                           grid_length=60,
                                           
                                           batch_size=9, #4, 6, 8 or 10
                                           
                                           theta_emulator=1/2,
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

# v90 <- pts9_6000$"All v90"
# Tc <- pts9_6000$"All Tc"
# values <- pts9_6000$"All values"
# hist(v90)
# hist(Tc)
# hist(values)
# plot(Tc,V90,
#      xlab="Tc",
#      ylab="V_{90}")

#dev.off()

#make a note of these!!!!!  
pts9_6000$"Best Design v90" #0.001169201
pts9_6000$"Best Design Tc" #5995.821
pts9_6000$"Best Design value" #-1.466055









# Figure 24 ---------------------------------------------------------------


#png(file="optimal_number_of_points.png",width=1200,height=1315,res=200)

plot(x=c(4,5,6,7,8,9),y=c(0.0008941561,
                          0.0006100015,
                          0.000475434,
                          0.0004916595,
                          0.0006746213,
                          0.001169201),
     xlab="Number of Points",
     ylab="V90",
     main="Finding the Optimal Number of Points",
     pch = 19
)
lines(x=c(4,5,6,7,8,9),y=c(0.0008941561,
                           0.0006100015,
                           0.000475434,
                           0.0004916595,
                           0.0006746213,
                           0.001169201))

#dev.off()






# End of Section 7 --------------------------------------------------------
