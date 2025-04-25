
# Section 2: Introduction to Emulation ------------------------------------

# Figures 1-9

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

##############################################################################################.
### Define simple Bayes Linear emulator for single input in 2D ###

simple_BL_emulator_v2 <- function(x,              # the emulator prediction point
                                  xD,             # the run input locations xD
                                  D,              # the run outputs D = (f(x^1),...,f(x^n))
                                  theta = 1,      # the correlation lengths
                                  sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])
                                  E_f = 0         # prior expectation of f: E(f(x)) = 0 
){
  
  # store length of runs D  
  n <- length(D) #number of runs we have from true computer model
  
  ### Define Covariance structure of f(x): Cov[f(x),f(xdash)] ###
  # Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-(x-xdash)^2/theta^2)    # XXX Old 1D version
  Cov_fx_fxdash <- function(x,xdash) sigma^2 * exp(-sum((x-xdash)^2)/theta^2) # XXX New 2D version
  #squares each component then sums them so give Euclidean distance
  
  ### Define 5 objects needed for BL adjustment ###
  # Create E[D] vector
  E_D <- rep(E_f,n)
  
  # Create Var_D matrix:
  Var_D <- matrix(0,nrow=n,ncol=n)
  # for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i],xD[j])  # XXX Old 1D version
  for(i in 1:n) for(j in 1:n) Var_D[i,j] <- Cov_fx_fxdash(xD[i,],xD[j,])  # XXX New 2D version
  
  # Create E[f(x)]
  E_fx <- E_f
  
  # Create Var_f(x) 
  Var_fx <- sigma^2
  
  # Create Cov_fx_D row vector
  Cov_fx_D <- matrix(0,nrow=1,ncol=n)
  # for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j])    # XXX Old 1D version
  for(j in 1:n) Cov_fx_D[1,j] <- Cov_fx_fxdash(x,xD[j,])    # XXX New 2D version
  #covariance between our new point and jth run input point
  #note that both x and xD[j,] are vectors of length 2
  
  ### Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x) ###
  ED_fx   <-  E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)   
  VarD_fx <-  Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)  
  #t(Cov_fx_D) = Cov[D,f(x)]
  
  ### return emulator expectation and variance ###
  return(c("ExpD_f(x)"=ED_fx,"VarD_f(x)"=VarD_fx))  
  
}
### End simple Bayes Linear emulator for single input in 2D ###
##############################################################################################.

### define filled contour plot function for emulator output ###
emul_fill_cont <- function(
    cont_mat,            # matrix of values we want contour plot of 
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






# Figure 1 ----------------------------------------------------------------
### Define actual computer model/simulator ###
f <- function(x) cos(3*pi*x) #function named f with one input x
#f <- function(x) sin(2*pi*x)

### Define run locations ###
xD <- seq(0,1,0.2) #[0,0.2,0.4,0.6,0.8,1]

### Perform 6 runs of model and store as D (this would take days for realistic 
# example!) ###
D <- f(xD)


xP <- seq(0.001,0.999,len=201)

### Evaluate emulator over 201 prediction points xP ###
em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))  # t(): gives transpose of the matrix


#png(file="Example_1D_Deterministic_Emulator.png",width=2400,height=1600,res=200)

### Plot the emulator
plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Emulator Output")

#dev.off() 





# Figure 2 ----------------------------------------------------------------
### Define actual 2D computer model/simulator ###
f <- function(x) 1.2*cos(2*pi*(x[,1]+0.3)) - 0.7*cos(2*pi*(x[,1]+0.3)*(x[,2]-0.2))

### Define run locations ###
D_grid <- c(0.05,0.35,0.65,0.95)
xD <- as.matrix(expand.grid("x1"=D_grid,"x2"=D_grid)) #these are our input pts
#                                                      to be evaluated by true fn
# 4×4 grid design 
#       {0.05,0.35,0.65,0.95}×{0.05,0.35,0.65,0.95}

# Now we can directly evaluate these runs using 𝑓(𝑥):
### Perform 16 runs of model and store as D (this would takes days for realistic example!) ###
D <- f(xD)

### Define 50x50 grid of prediction points xP for emulator evaluation ###
x_grid <- seq(-0.001,1.001,len=50)
xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))

### Evaluate emulator over 50x50=2500 prediction points xP ###
em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=0.45,sigma=1,E_f=0))   

### store emulator output as matrices to aid plotting ###
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

library(viridisLite)
exp_cols <- magma #see ?viridis
var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
diag_cols <- turbo #see ?viridis

#png(file="BIGGER_Example_2D_Deterministic_Emulator_Expectation.png",width=1400,height=1350,res=200)

emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,        # this sets the colour scheme
               main="Emulator Adjusted Expectation E_D[f(x)]")

#dev.off()

#png(file="BIGGER_Example_2D_Deterministic_Emulator_Variance.png",width=1400,height=1350,res=200)

emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
               color.palette=var_cols,
               main="Emulator Adjusted Variance Var_D[f(x)]")

#dev.off()

### Evaluate true function and store in matrix for diagnostic comparisons ###
fxP_mat <- matrix(f(xP),nrow=length(x_grid),ncol=length(x_grid)) 

#png(file="BIGGER_Example_2D_Deterministic_Emulator_True.png",width=1400,height=1350,res=200)

emul_fill_cont(cont_mat=fxP_mat,cont_levs=seq(-2,2,0.2),xD=xD,x_grid=x_grid,
               color.palette=exp_cols,
               main="True Computer Model Function f(x)")

#dev.off()

### Evaluate diagnostics S_D(x) and store in matrix ###
S_diag_mat <- (E_D_fx_mat - fxP_mat) / sqrt(Var_D_fx_mat)

#png(file="BIGGER_Example_2D_Deterministic_Emulator_Diagnostics.png",width=1400,height=1350,res=200)

emul_fill_cont(cont_mat=S_diag_mat,cont_levs=seq(-3,3,0.25),xD=xD,x_grid=x_grid,
               xD_col="purple",
               color.palette=diag_cols,
               main="Emulator Diagnostics S_D[f(x)]")

#dev.off()




# Figure 3 ----------------------------------------------------------------
### Define actual computer model/simulator ###
f <- function(x) cos(3*pi*x)

### Define run locations ###
xD <- seq(0,1,0.2) #[0,0.2,0.4,0.6,0.8,1]

### Perform 6 runs of model and store as D (this would take days for realistic 
# example!) ###
D <- f(xD)

xP <- seq(0.001,0.999,len=201)

#putting various prior expectations on same graph
##############################################################################################.
### Function to plot simple emulator output
allonsamegraph_plot_BL_emulator_V1 <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))^T
    maintitle=NULL  # title of the plot, blank if NULL
){
  
  ### plot emulator output ###
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.2,1.2),ty="l",col="yellow",lwd=2.5,lty=2,
       xlab="Input parameter x",ylab="Output f(x)",main=maintitle) # main argument: plot title (INPUT TO SPECIFY)
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="yellow",lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="yellow",lwd=2.5)
  
  ### plot true function: we would not normally be able to do this! ###
  lines(xP,f(xP),lwd=2,lty=1)
  
  ### Plot the runs ### 
  points(xD,D,pch=21,col=1,bg="green",cex=1.5) #PLOT THE TRUE RUNS OF THE EXPENSIVE MODEL
  legend('topright',legend=c("True function f(x)",
                             "Model Evaluations",
                             "Emulator Expectation: E[f(x)] = -2",
                             "Emulator Expectation: E[f(x)] = -1",
                             "Emulator Expectation: E[f(x)] = 0",
                             "Emulator Expectation: E[f(x)] = 1",
                             "Emulator Expectation: E[f(x)] = 2"),
         lty=c(1,NA,2,2,2,2,2),pch=c(NA,16,NA,NA,NA,NA,NA),col=c(1,"green","yellow","orange","deeppink","purple","cadetblue1"),lwd=2.5,pt.cex=1.3)
}
### Function to plot simple emulator output
##############################################################################################.

#yellow orange deeppink purple


em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=-2))

prior_expectation_seq <- c(-1,0,1,2)
colours <- c("orange","deeppink","purple","cadetblue1")

#png(file="Varying_Prior_Expectation_allononegraph.png",width=2800,height=1600,res=200)

### Plot emulator output in each case, note use of "paste" for plot title ###
allonsamegraph_plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle=paste("Emulator Output Varying Prior Expectation"))

### for loop over different sigma values in sigma_seq ###
for(i in 1:length(prior_expectation_seq)){
  ### Evaluate emulator over 201 prediction points xP with E_f=prior_expectation_seq[i] ###
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=prior_expectation_seq[i]))
  
  lines(xP,em_out[,"ExpD_f(x)"],col=colours[i],lwd=2.5,lty=2)
  
  ### Plot emulator prediction intervals in each case, note use of "paste" for plot title ###
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col=colours[i],lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col=colours[i],lwd=2.5)
}

legend('bottomleft',legend=c("Pred. Interval: E[f(x)] = -2",
                             "Pred. Interval: E[f(x)] = -1",
                             "Pred. Interval: E[f(x)] = 0",
                             "Pred. Interval: E[f(x)] = 1",
                             "Pred. Interval: E[f(x)] = 2"),
       lty=c(1,1,1,1,1),col=c("yellow","orange","deeppink","purple","cadetblue1"),lwd=2.5)

lines(xP,f(xP),lwd=2,lty=1)

### Plot the runs ### 
points(xD,D,pch=21,col=1,bg="green",cex=1.5) #PLOT THE TRUE RUNS OF THE EXPENSIVE MODEL


#dev.off()






# Figure 4 ----------------------------------------------------------------
### Define actual computer model/simulator ###
f <- function(x) cos(3*pi*x)

### Define run locations ###
xD <- seq(0,1,0.2) #[0,0.2,0.4,0.6,0.8,1]

### Perform 6 runs of model and store as D (this would take days for realistic 
# example!) ###
D <- f(xD)

xP <- seq(0.001,0.999,len=201)

##############################################################################################.
### Function to plot simple emulator output
varying_prior_var_plot_BL_emulator_V1 <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))^T
    maintitle=NULL  # title of the plot, blank if NULL
){
  
  ### plot emulator output ###
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.3,1.3),ty="l",col="blue",lwd=2.5,
       xlab="Input parameter x",ylab="Output f(x)",main=maintitle) # main argument: plot title (INPUT TO SPECIFY)
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="yellow",lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="yellow",lwd=2.5)
  
  ### plot true function: we would not normally be able to do this! ###
  lines(xP,f(xP),lwd=2,lty=1)
  
  legend('topright',legend=c("Emulator Expectation",
                             "True function f(x)",
                             "Model Evaluations"),
         lty=c(1,1,NA),pch=c(NA,NA,16),col=c("blue",1,"green"),lwd=2.5,pt.cex=1.3)
}

### Function to plot simple emulator output
##############################################################################################.


# Now let’s vary the prior variance parameter 𝜎^2 = Var[𝑓(𝑥)]. 
# We will choose a series of values for 𝜎 and make a fresh emulator plot for 
# each choice:

### make sequence of sigma values to use for emulation ###
sigma_seq <- c(0.5,0.75,1)
colours <- c("orange","deeppink","purple")


#png(file="Varying_Prior_Variance.png",width=2400,height=1600,res=200)

em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.25,E_f=0))
varying_prior_var_plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Emulator Output Varying Sigma")

### for loop over different sigma values in sigma_seq ###
for(i in 1:length(sigma_seq)){
  ### Evaluate emulator over 201 prediction points xP with sigma=sigma_seq[i] ###
  em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=sigma_seq[i],E_f=0))
  
  ### Plot emulator prediction intervals in each case, note use of "paste" for plot title ###
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col=colours[i],lwd=2.5)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col=colours[i],lwd=2.5)
}

legend('bottomleft',legend=c("Pred. Interval: sigma = 0.25",
                             "Pred. Interval: sigma = 0.5",
                             "Pred. Interval: sigma = 0.75",
                             "Pred. Interval: sigma = 1"),
       lty=c(1,1,1,1),col=c("yellow","orange","deeppink","purple"),lwd=2.5)

### Plot the runs ### 
points(xD,D,pch=21,col=1,bg="green",cex=1.5) #PLOT THE TRUE RUNS OF THE EXPENSIVE MODEL

#dev.off()






# Figure 5 ----------------------------------------------------------------
### Define actual computer model/simulator ###
f <- function(x) cos(3*pi*x)

### Define run locations ###
xD <- seq(0,1,0.2) #[0,0.2,0.4,0.6,0.8,1]

### Perform 6 runs of model and store as D (this would take days for realistic 
# example!) ###
D <- f(xD)

xP <- seq(0.001,0.999,len=201)

#need nugget emulator as singular matrix
simple_BL_emulator_v1_with_nugget <- function(x,              # the emulator prediction point
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



##############################################################################################.
### Function to plot simple emulator output
corrlength_plot_BL_emulator_V1 <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))^T
    maintitle=NULL  # title of the plot, blank if NULL
){
  
  ### plot emulator output ###
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.8,1.8),ty="l",col="blue",lwd=2.5,
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

#1D example

# Now let’s vary the correlation length parameter 𝜃. 
# We will choose a series of values for 𝜃 and make a fresh emulator plot for 
# each choice:

### make sequence of theta values to use for emulation ###
theta_seq <- c(0.005,0.1,0.2,0.3,0.4,0.7)

#png(file="BIGGER_varyingtheta%d.png",width=2100,height=1400,res=200)

### for loop over different theta values in theta_seq ###
for(i in 1:length(theta_seq)){
  ### Evaluate emulator over 201 prediction points xP with theta=theta_seq[i] ###
  em_out <- t(sapply(xP,simple_BL_emulator_v1_with_nugget,xD=xD,D=D,theta=theta_seq[i],sigma=0.5,E_f=0))
  
  ### Plot emulator output in each case, note use of "paste" for plot title ###
  corrlength_plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle=paste("Theta =",theta_seq[i]))
}

#dev.off()






# Figure 6 ----------------------------------------------------------------
### Define actual 2D computer model/simulator ###
f <- function(x) 1.2*cos(2*pi*(x[,1]+0.3)) - 0.7*cos(2*pi*(x[,1]+0.3)*(x[,2]-0.2))

### Define run locations ###
D_grid <- c(0.05,0.35,0.65,0.95)
xD <- as.matrix(expand.grid("x1"=D_grid,"x2"=D_grid)) #these are our input pts
#                                                      to be evaluated by true fn
# 4×4 grid design 
#       {0.05,0.35,0.65,0.95}×{0.05,0.35,0.65,0.95}

# Now we can directly evaluate these runs using 𝑓(𝑥):
### Perform 16 runs of model and store as D (this would takes days for realistic example!) ###
D <- f(xD)

### Define 50x50 grid of prediction points xP for emulator evaluation ###
x_grid <- seq(-0.001,1.001,len=50)
xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))

library(viridisLite)
exp_cols <- magma #see ?viridis
var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
diag_cols <- turbo #see ?viridis

#2D example

### Define a vector of theta values to use ###
theta_seq <- c(0.05,0.1,0.15,0.2,0.35,0.4,0.45)

#png(file="BIGGER_2Dvaryingtheta%d.png",width=1400,height=1350,res=200)

### loop over vector of theta values ###
for(i in 1:length(theta_seq)){
  
  ### Evaluate emulator over 201 prediction points xP and store in matrices ###
  em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=theta_seq[i],sigma=1,E_f=0))   
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  ### Plot filled contour plot of emulator expectation ###
  emul_fill_cont(cont_mat=Var_D_fx_mat,nlev=12,xD=xD,x_grid=x_grid,color.palette=var_cols,
                 main=paste("Emulator Adjusted Variance Var_D[f(x)], theta =",theta_seq[i]))
}

#dev.off()







# Figure 7 ----------------------------------------------------------------
### Define actual 2D computer model/simulator ###
f <- function(x) 1.2*cos(2*pi*(x[,1]+0.3)) - 0.7*cos(2*pi*(x[,1]+0.3)*(x[,2]-0.2))

### Define run locations ###
D_grid <- c(0.05,0.35,0.65,0.95)
xD <- as.matrix(expand.grid("x1"=D_grid,"x2"=D_grid)) #these are our input pts
#                                                      to be evaluated by true fn
# 4×4 grid design 
#       {0.05,0.35,0.65,0.95}×{0.05,0.35,0.65,0.95}

# Now we can directly evaluate these runs using 𝑓(𝑥):
### Perform 16 runs of model and store as D (this would takes days for realistic example!) ###
D <- f(xD)

### Define 50x50 grid of prediction points xP for emulator evaluation ###
x_grid <- seq(-0.001,1.001,len=50)
xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))

library(viridisLite)
exp_cols <- magma #see ?viridis
var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
diag_cols <- turbo #see ?viridis

#2D example to show regions must overlap


### Define a vector of theta values to use ###
new_theta_seq <- c(0.1,0.15,0.35)

#png(file="2nd_BIGGER_2Dvaryingtheta_so_regions_overlap%d.png",width=1400,height=1350,res=200)

par(mar = c(5.1, 3.1, 4.1, 2.1))

### loop over vector of theta values ###
for(i in 1:length(new_theta_seq)){
  
  ### Evaluate emulator over 201 prediction points xP and store in matrices ###
  em_out <- t(apply(xP,1,simple_BL_emulator_v2,xD=xD,D=D,theta=new_theta_seq[i],sigma=1,E_f=0))   
  
  #extracting the emulator expectation and variance
  E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  
  
  #graphing
  emul_fill_cont(cont_mat=Var_D_fx_mat,nlev=12,xD=xD,x_grid=x_grid,color.palette=var_cols,
                 main=paste("Emulator Adjusted Variance Var_D[f(x)], theta =",new_theta_seq[i]))
  emul_fill_cont(cont_mat=E_D_fx_mat,cont_levs=NULL,xD=xD,x_grid=x_grid,
                 color.palette=exp_cols,main=paste("Emulator Adjusted Expectation E_D[f(x)], theta =",new_theta_seq[i]))
}

#dev.off()





# Figure 8 ----------------------------------------------------------------
### Define actual computer model/simulator ###
f <- function(x) cos(3*pi*x)

### Define run locations ###
xD <- seq(0,1,0.2) #[0,0.2,0.4,0.6,0.8,1]

### Perform 6 runs of model and store as D (this would take days for realistic 
# example!) ###
D <- f(xD)

xP <- seq(0.001,0.999,len=201)

#Large nugget - stochasticity
#one with nugget large 
#another with nugget very large

nugget_simple_BL_emulator_v1 <- function(x,              # the emulator prediction point
                                         xD,             # the run input locations xD
                                         D,              # the run outputs D = (f(x^1),...,f(x^n))
                                         theta = 1,      # the correlation lengths     (1 if not specified otherwise)
                                         sigma = 1,      # the prior SD sigma sqrt(Var[f(x)])     (1 if not specified otherwise)
                                         E_f = 0,        # prior expectation of f: E(f(x)) = 0     (0 if not specified otherwise)
                                         nugget = 10^(-6)
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
  
  #nugget part
  diagonal <- diag(nugget,n,n)
  Var_D <- Var_D + sigma^2*diagonal
  #end of nugget part
  
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

stochastic_plot_BL_emulator_V1 <- function(
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
  
  
  ### Plot the runs ### 
  points(xD,D,pch=21,col=1,bg="green",cex=1.5) #PLOT THE TRUE RUNS OF THE EXPENSIVE MODEL
  legend('topright',legend=c("Emulator Expectation","Emulator Prediction Interval","Model Evaluations"),
         lty=c(1,1,NA),pch=c(NA,NA,16),col=c("blue","red","green"),lwd=2.5,pt.cex=1.3)
}                                 

em_out <- t(sapply(xP,nugget_simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0, nugget=0.05^2))  
em_out2 <- t(sapply(xP,nugget_simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0, nugget=0.2^2))  


#png(file="BIGGER_large_nugget_%d.png",width=2100,height=1400,res=200)

### Plot the emulator
stochastic_plot_BL_emulator_V1(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Emulator Output, Nugget = 0.05")
stochastic_plot_BL_emulator_V1(em_out=em_out2,xP=xP,xD=xD,D=D,maintitle="Emulator Output, Nugget = 0.2")

#dev.off()





# Figure 9 ----------------------------------------------------------------
### Define actual computer model/simulator ###
f <- function(x) cos(3*pi*x)

### Define run locations ###
xD <- seq(0,1,0.2) #[0,0.2,0.4,0.6,0.8,1]

### Perform 6 runs of model and store as D (this would take days for realistic 
# example!) ###
D <- f(xD)

xP <- seq(0.001,0.999,len=201)

#Graph showing poor at extrapolation
#Does not extrapolate well outside range of runs. Resorts back to prior.

f <- function(x) cos(3*pi*x)
xD <- seq(0,1,0.2) #same runs
D <- f(xD)
xP <- seq(-0.999,1.999,len=601) #new range from -1 to 2


### Function to plot simple emulator output
plot_BL_emulator_V1_BIGGER <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))^T
    maintitle=NULL  # title of the plot, blank if NULL
){
  
  ### plot emulator output ###
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.9,2.4),ty="l",col="blue",lwd=2.5,
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

#png(file="BIGGER_badatextrapolation.png",width=1950,height=1625,res=200)

em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=0))
plot_BL_emulator_V1_BIGGER(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Prior Expectation = 0")

#dev.off()




expectation2_plot_BL_emulator_V1_BIGGER <- function(
    em_out,         # a nP x 2 matrix of emulator outputs
    xP,             # vector of nP inputs where emulator evaluated
    xD,             # the run input locations
    D,              # the run outputs D = (f(x^1),...,f(x^n))^T
    maintitle=NULL  # title of the plot, blank if NULL
){
  
  ### plot emulator output ###
  plot(xP,em_out[,"ExpD_f(x)"],ylim=c(-1.1,4.4),ty="l",col="blue",lwd=2.5,
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

#png(file="BIGGER_badatextrapolation_expectation2.png",width=1950,height=1625,res=200)
em_out <- t(sapply(xP,simple_BL_emulator_v1,xD=xD,D=D,theta=0.25,sigma=0.5,E_f=2))
expectation2_plot_BL_emulator_V1_BIGGER(em_out=em_out,xP=xP,xD=xD,D=D,maintitle="Prior Expectation = 2")
#dev.off()





# End of Section 2 --------------------------------------------------------