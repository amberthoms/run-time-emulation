
# Section 10: Second Wave of the Iterative History Match ------------------

# Figures 32-34

# Functions ---------------------------------------------------------------



# Other: loading the new evaluations etc. -------------------------------------------------------------------

## Stochasticity Runs ------------------------------------------------------

stochasticity_data <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stochasticity_tests/Runtimes.txt")
log10_stochasticity_data <- log10(stochasticity_data$V1)

## New Data ----------------------------------------------------------------

new_runs <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/LATEST FOLDER/OneDrive_1_04-04-2025/Stellar_AGN_feedback_2D_v2/simulation_overview.txt")
matrix_new_runs <- as.matrix(new_runs)
colnames(matrix_new_runs) <- c("Emulator Coordinate x", 
                               "Emulator Coordinate y", 
                               "Runtime [CPU-Hours]", 
                               "Batch Number", 
                               "fmin",	
                               "fmax",	
                               "rho0 [cm^-3]",
                               "T_AGN [K]")

# New runs based on our best configuration for run time
x_old_new_runs <- head(matrix_new_runs[,"Emulator Coordinate x"],-1)
y_old_new_runs <- head(matrix_new_runs[,"Emulator Coordinate y"],-1)
xD_old_new_runs <- matrix(data=c(x_old_new_runs,y_old_new_runs),ncol=2,byrow = FALSE)
D_old_new_runs <- log10(head(matrix_new_runs[,"Runtime [CPU-Hours]"],-1))
# We will use mean(D_old_new_runs) as our E_f_val
# We will use sd(D_old_new_runs) for our sigma
# Recall we calculate our nugget using:
sd(log10_stochasticity_data)/sd(D_old_new_runs) # 0.1203885
# therefore use 0.1203885^2 as our nugget term

# With the additional odd run in the bottom left corner
x_old_new_extra_runs <- matrix_new_runs[,"Emulator Coordinate x"]
y_old_new_extra_runs <- matrix_new_runs[,"Emulator Coordinate y"]
xD_old_new_extra_runs <- matrix(data=c(x_old_new_extra_runs,y_old_new_extra_runs),ncol=2,byrow = FALSE)
D_old_new_extra_runs <- log10(matrix_new_runs[,"Runtime [CPU-Hours]"])
# We will use mean(D_old_new_extra_runs) as our E_f_val
# We will use sd(D_old_new_extra_runs) for our sigma
# Recall we calculate our nugget using:
sd(log10_stochasticity_data)/sd(D_old_new_extra_runs) # 0.1176641
# therefore use 0.1176641^2 as our nugget term






# Figure 32 ---------------------------------------------------------------


## 32a ---------------------------------------------------------------------

# See Figure 27 in Section 8.3.2


## 32b ---------------------------------------------------------------------


# See Figure 30 in Section 9.1 for function for simulating total run times

## Add vertical line showing true run time on histogram of simulated run times --------

new_runs_indices <- which(matrix_new_runs[,"Batch Number"]==1)
actual_time <- sum(matrix_new_runs[,"Runtime [CPU-Hours]"][new_runs_indices])
actual_time #4167.24

winning_design_simulations <- readRDS("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/My work/RStudio/week 19/winningdesignhist1_results.RData")

#png(file="actual_time_line_simulated_histogram_final.png",width=2000,height=1200,res=200)
hist(winning_design_simulations$`simulated totals`,
     xlab="Total Run Time (natural scale)",
     main = "Histogram of 1000 Simulated Total Run Times")
abline(v=c(winning_design_simulations$`expected total`,winning_design_simulations$`mean total`,quantile(winning_design_simulations$`simulated totals`,c(0.025,0.975)),actual_time),
       col=c("red","gold","green","blue","hotpink"),
       lty=c(1,2,1,1,2),
       lwd=2.5)
legend(x="topright",
       legend=c("Expectation","Mean","2.5th Percentile","97.5th Percentile","Actual Run Time"),
       lty=c(1,2,1,1,2),
       col=c("red","gold","green","blue","hotpink"),
       lwd=2.5)
#dev.off()





# Figure 33 ---------------------------------------------------------------


## History Matching --------------------------------------------------------


## functions required ------------------------------------------------------

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
  if(using_pdist)  Cov_fx_D <- Cov_fx_fxdash( as.matrix(pdist::pdist(xP,xD)) )   # XXX New V3
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

### define filled contour plot function for emulator output ###
emul_fill_cont_V2 <- function(
    cont_mat,            # matrix of values we want contour plot of 
    cont_levs=NULL,      # contour levels (NULL: automatic selection)
    cont_levs_lines=NULL,   # contour levels for lines (NULL: automatic selection)
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
                   if(is.null(cont_levs_lines)) contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs,lwd=0.8) # plot usual contour lines 
                   if(!is.null(cont_levs_lines)) {
                     contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs,lwd=0.4,labels="")   # plot thin contour lines 
                     contour(x_grid,x_grid,cont_mat,add=TRUE,levels=cont_levs_lines,lwd=2)   # plot thick contour lines
                   }
                   if(plot_xD) points(xD,pch=21,col=1,bg=xD_col,cex=1.5)})  # plot design points
}


exp_cols <- magma
var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
diag_cols <- turbo


hm_plot <- function(hm_z,
                    hm_sigma_e,
                    hm_sigma_epsilon,
                    hm_grid_length,
                    hm_old_vec,
                    hm_new_vec,
                    hm_D_runs,
                    hm_theta,
                    hm_sigma,
                    hm_E_f,
                    hm_nugget
){
  z <- hm_z
  sigma_e <- hm_sigma_e
  sigma_epsilon <- hm_sigma_epsilon
  x_grid <- seq(-0.001,1.001,len=hm_grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  xD <- matrix(hm_old_vec,ncol=2,byrow = TRUE)
  D <- hm_D_runs
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=hm_theta,sigma=hm_sigma,E_f=hm_E_f,nugget=hm_nugget)  
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Imp_mat <- sqrt( (E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2) ) 
  
  imp_cols <- function(n) turbo(n,begin=0.15,end=1)
  imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),20,40,61)
  emul_fill_cont_V2(cont_mat=Imp_mat,
                    cont_levs=imp_levs,
                    cont_levs_lines=3,
                    xD=xD,
                    x_grid=x_grid,
                    xD_col="purple",
                    color.palette=imp_cols,main="Implausibility I(x)")
  
}



## new variables required for hm -------------------------------------------

z_10.5 <- 0.00512
sigma_e_10.5 <- 0.000512
sigma_epsilon_10.5 <- 0.000512

## Original Latin Hypercube ------------------------------------------------

D_10.5 <- c(0.0145267104,
            0.0012496095,
            0.0000000000,
            0.0174945330,
            0.0053108404,
            0.0006248047,
            0.0153077163,
            0.0110902843,
            0.0010934083)
sd_10.5 <- sd(D_10.5)
mean_10.5 <- mean(D_10.5)
nugget_10.5 <- 0.03088729^2

#png(file="batch0_hm_%d.png",width=1400,height=1200,res=200)
hm_plot(hm_z=z_10.5,
        hm_sigma_e=sigma_e_10.5,
        hm_sigma_epsilon=sigma_epsilon_10.5,
        hm_grid_length=60,
        hm_old_vec=c(t(xD)),
        hm_D_runs=D_10.5,
        hm_theta=1/2,
        hm_sigma=sd_10.5,
        hm_E_f=mean_10.5,
        hm_nugget=nugget_10.5)
#dev.off()


## ASIDE -------------------------------------------------------------------

# Before we can do hm for new batches we need to actually calculate the 10.5
# stellar mass function at each of the new evaluations!!!!!!!

Run_9 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/LATEST FOLDER/OneDrive_1_04-04-2025/Stellar_AGN_feedback_2D_v2/Run_9/halo_catalogue.txt")
Run_10 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/LATEST FOLDER/OneDrive_1_04-04-2025/Stellar_AGN_feedback_2D_v2/Run_10/halo_catalogue.txt")
Run_11 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/LATEST FOLDER/OneDrive_1_04-04-2025/Stellar_AGN_feedback_2D_v2/Run_11/halo_catalogue.txt")
Run_12 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/LATEST FOLDER/OneDrive_1_04-04-2025/Stellar_AGN_feedback_2D_v2/Run_12/halo_catalogue.txt")
Run_13 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/LATEST FOLDER/OneDrive_1_04-04-2025/Stellar_AGN_feedback_2D_v2/Run_13/halo_catalogue.txt")
Run_14 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/LATEST FOLDER/OneDrive_1_04-04-2025/Stellar_AGN_feedback_2D_v2/Run_14/halo_catalogue.txt")
Run_15 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/LATEST FOLDER/OneDrive_1_04-04-2025/Stellar_AGN_feedback_2D_v2/Run_15/halo_catalogue.txt")

#we add 10 to log10(...) which is equivalent to log10(...x10^10)
hist_run_9 <- hist(log10(Run_9$V2)+10, 
                   plot=F,
                   breaks=seq(6.5,12,by=0.5))

hist_run_10 <- hist(log10(Run_10$V2)+10, 
                    plot=F,
                    breaks=seq(6.5,12,by=0.5))

hist_run_11 <- hist(log10(Run_11$V2)+10, 
                    plot=F,
                    breaks=seq(6.5,12,by=0.5))

hist_run_12 <- hist(log10(Run_12$V2)+10, 
                    plot=F,
                    breaks=seq(6.5,12,by=0.5))

hist_run_13 <- hist(log10(Run_13$V2)+10, 
                    plot=F,
                    breaks=seq(6.5,12,by=0.5))

hist_run_14 <- hist(log10(Run_14$V2)+10, 
                    plot=F,
                    breaks=seq(6.5,12,by=0.5))

hist_run_15 <- hist(log10(Run_15$V2)+10, 
                    plot=F,
                    breaks=seq(6.5,12,by=0.5))



#number of galaxies per log_10 (M_st)
run9_step1 <- hist_run_9$counts/0.5
run10_step1 <- hist_run_10$counts/0.5
run11_step1 <- hist_run_11$counts/0.5
run12_step1 <- hist_run_12$counts/0.5
run13_step1 <- hist_run_13$counts/0.5
run14_step1 <- hist_run_14$counts/0.5
run15_step1 <- hist_run_15$counts/0.5


#number of galaxies per log10(M_st) per unit volume
run9_step2 <- run9_step1/6402
run10_step2 <- run10_step1/6402
run11_step2 <- run11_step1/6402
run12_step2 <- run12_step1/6402
run13_step2 <- run13_step1/6402
run14_step2 <- run14_step1/6402
run15_step2 <- run15_step1/6402


# extract y values at x=10.5
library(pracma)

run_9_10.5 <- interp1(x=hist_run_9$mids,
                      y=run9_step2,
                      xi=10.5,
                      method="linear")
run_10_10.5 <- interp1(x=hist_run_9$mids,
                       y=run10_step2,
                       xi=10.5,
                       method="linear")
run_11_10.5 <- interp1(x=hist_run_9$mids,
                       y=run11_step2,
                       xi=10.5,
                       method="linear")
run_12_10.5 <- interp1(x=hist_run_9$mids,
                       y=run12_step2,
                       xi=10.5,
                       method="linear")
run_13_10.5 <- interp1(x=hist_run_9$mids,
                       y=run13_step2,
                       xi=10.5,
                       method="linear")
run_14_10.5 <- interp1(x=hist_run_9$mids,
                       y=run14_step2,
                       xi=10.5,
                       method="linear")
run_15_10.5 <- interp1(x=hist_run_9$mids,
                       y=run15_step2,
                       xi=10.5,
                       method="linear")


extracted_at_10.5 <- c(run_9_10.5,
                       run_10_10.5,
                       run_11_10.5,
                       run_12_10.5,
                       run_13_10.5,
                       run_14_10.5,
                       run_15_10.5)


# write out so easier (and do not need to run above every time)
D_batch1_minus1 <- c(0.0059356451,
                     0.0000000000,
                     0.0000000000,
                     0.0000000000,
                     0.0006248047,
                     0.0103092784,
                     0.0189003436)

all_D <- c(D_10.5,D_batch1_minus1)
all_except_extra <- head(all_D,-1)

# and then calculate the nuggets!!!

# extract stochastic y values at x=10.5
# see observeddata and week16_forweek17_hm
stochastic_extracted_at_10.5 <- c(0.005154639,
                                  0.005467042,
                                  0.004998438,
                                  0.005154639,
                                  0.005467042,
                                  0.005154639,
                                  0.005467042,
                                  0.004998438,
                                  0.004842237,
                                  0.005310840)


#nugget calculations:
sd(stochastic_extracted_at_10.5)/sd(all_except_extra) # 0.03435351

sd(stochastic_extracted_at_10.5)/sd(all_D) # 0.03135506


## New runs based on our best configuration for run time -------------------

#new hmplot so wave 2 in hot pink
hm_plot <- function(hm_z,
                    hm_sigma_e,
                    hm_sigma_epsilon,
                    hm_grid_length,
                    hm_old_vec,
                    hm_new_vec,
                    hm_D_runs,
                    hm_theta,
                    hm_sigma,
                    hm_E_f,
                    hm_nugget
){
  z <- hm_z
  sigma_e <- hm_sigma_e
  sigma_epsilon <- hm_sigma_epsilon
  x_grid <- seq(-0.001,1.001,len=hm_grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  xD <- matrix(hm_old_vec,ncol=2,byrow = TRUE)
  D <- hm_D_runs
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=hm_theta,sigma=hm_sigma,E_f=hm_E_f,nugget=hm_nugget)  
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Imp_mat <- sqrt( (E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2) ) 
  
  imp_cols <- function(n) turbo(n,begin=0.15,end=1)
  imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),20,40,61)
  emul_fill_cont_V2(cont_mat=Imp_mat,
                    cont_levs=imp_levs,
                    cont_levs_lines=3,
                    xD=xD,
                    x_grid=x_grid,
                    xD_col=c(rep("purple",9),rep("hotpink",6)),
                    color.palette=imp_cols,main="Implausibility I(x)")
  
}


D_10.5_old_new_runs <- all_except_extra
sd_10.5_old_new_runs <- sd(D_10.5_old_new_runs)
mean_10.5_old_new_runs <- mean(D_10.5_old_new_runs)
nugget_10.5_old_new_runs <- 0.03435351^2

#png(file="different_cols_batch0_and_1_hm_%d.png",width=1400,height=1200,res=200)
hm_plot(hm_z=z_10.5,
        hm_sigma_e=sigma_e_10.5,
        hm_sigma_epsilon=sigma_epsilon_10.5,
        hm_grid_length=60,
        hm_old_vec=c(t(xD_old_new_runs)),
        hm_D_runs=D_10.5_old_new_runs,
        hm_theta=1/2,
        hm_sigma=sd_10.5_old_new_runs,
        hm_E_f=mean_10.5_old_new_runs,
        hm_nugget=nugget_10.5_old_new_runs)
#dev.off()



## With the additional odd run in the bottom left corner -------------------

#new hmplot so wave 2 in hot pink and wave 3 in blue
hm_plot <- function(hm_z,
                    hm_sigma_e,
                    hm_sigma_epsilon,
                    hm_grid_length,
                    hm_old_vec,
                    hm_new_vec,
                    hm_D_runs,
                    hm_theta,
                    hm_sigma,
                    hm_E_f,
                    hm_nugget
){
  z <- hm_z
  sigma_e <- hm_sigma_e
  sigma_epsilon <- hm_sigma_epsilon
  x_grid <- seq(-0.001,1.001,len=hm_grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  xD <- matrix(hm_old_vec,ncol=2,byrow = TRUE)
  D <- hm_D_runs
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=hm_theta,sigma=hm_sigma,E_f=hm_E_f,nugget=hm_nugget)  
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Imp_mat <- sqrt( (E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2) ) 
  
  imp_cols <- function(n) turbo(n,begin=0.15,end=1)
  imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),20,40,61)
  emul_fill_cont_V2(cont_mat=Imp_mat,
                    cont_levs=imp_levs,
                    cont_levs_lines=3,
                    xD=xD,
                    x_grid=x_grid,
                    xD_col=c(rep("purple",9),rep("hotpink",6),"cyan"),
                    color.palette=imp_cols,main="Implausibility I(x)")
  
}



D_10.5_old_new_extra_runs <- all_D
sd_10.5_old_new_extra_runs <- sd(D_10.5_old_new_extra_runs)
mean_10.5_old_new_extra_runs <- mean(D_10.5_old_new_extra_runs)
nugget_10.5_old_new_extra_runs <- 0.03135506^2

#png(file="different_cols_batch0_and_1_and_extra_hm_%d.png",width=1400,height=1200,res=200)
hm_plot(hm_z=z_10.5,
        hm_sigma_e=sigma_e_10.5,
        hm_sigma_epsilon=sigma_epsilon_10.5,
        hm_grid_length=60,
        hm_old_vec=c(t(xD_old_new_extra_runs)),
        hm_D_runs=D_10.5_old_new_extra_runs,
        hm_theta=1/2,
        hm_sigma=sd_10.5_old_new_extra_runs,
        hm_E_f=mean_10.5_old_new_extra_runs,
        hm_nugget=nugget_10.5_old_new_extra_runs)
#dev.off()












# Figure 34 ---------------------------------------------------------------


## 34a ---------------------------------------------------------------------

# See Figure 22a in Section 7.2

## 34b ---------------------------------------------------------------------

## Remake run time emulator ------------------------------------------------


## functions required ------------------------------------------------------

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
  if(using_pdist)  Cov_fx_D <- Cov_fx_fxdash( as.matrix(pdist::pdist(xP,xD)) )   # XXX New V3
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


## Original Latin Hypercube ------------------------------------------------
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


## New runs based on our best configuration for run time -------------------

#png(file="old_new_runs_emulator_log_runtime_%d.png",width=1400,height=1200,res=200)
nugget_argument_efficient_nugget_all_graphs(xD_matrix = xD_old_new_runs,
                                            D_vector = D_old_new_runs,
                                            grid_length = 50,
                                            theta_val = 1, #pick a large theta BECAUSE WE EXPECT RUN TIME FUNCTION TO BE SMOOTH!
                                            sigma_val = sd(D_old_new_runs), 
                                            E_f_val = mean(D_old_new_runs),
                                            
                                            nugget_val = 0.1203885^2) #must input it as squared as dealing with sigma squared
#dev.off()


## With the additional odd run in the bottom left corner -------------------

#png(file="old_new_extra_runs_plus_extra_emulator_log_runtime_%d.png",width=1400,height=1200,res=200)
nugget_argument_efficient_nugget_all_graphs(xD_matrix = xD_old_new_extra_runs,
                                            D_vector = D_old_new_extra_runs,
                                            grid_length = 50,
                                            theta_val = 1, #pick a large theta BECAUSE WE EXPECT RUN TIME FUNCTION TO BE SMOOTH!
                                            sigma_val = sd(D_old_new_extra_runs), 
                                            E_f_val = mean(D_old_new_extra_runs),
                                            
                                            nugget_val = 0.1176641^2) #must input it as squared as dealing with sigma squared
#dev.off()






# End of Section 10 -------------------------------------------------------
