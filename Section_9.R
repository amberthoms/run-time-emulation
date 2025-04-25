
# Section 9: Run Time Uncertainty -----------------------------------------

# Figures 29-31

# Functions ---------------------------------------------------------------

# Figure 29 ---------------------------------------------------------------


## First run this ----------------------------------------------------------

rm(list = ls())

#need some new small functions first

V90W2_function <- function(grid_length_val3,
                           
                           oldpts_vec_val3,
                           newpts_vec_val3,
                           
                           theta_val3,
                           sigma_val3,
                           E_f_val3,
                           nugget_val3,
                           
                           z_val3,
                           sigma_e_val3,
                           sigma_epsilon_val3,
                           
                           percentile_val3,
                           
                           val3_actual_stellarmassfunctions){
  
  # xD_here <- xD_val3
  # D_here <- val3_actual_stellarmassfunctions
  # x_grid <- seq(-0.001,1.001,len=grid_length_val3)
  # xP_here <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  # 
  # em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP_here,xD=xD_here,D=D_here,theta=theta_val3,sigma=sigma_val3,E_f=E_f_val3,nugget=nugget_val3)  
  # E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  # Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
  # 
  # Imp_mat <- sqrt( (E_D_fx_mat - z_val3)^2 / (Var_D_fx_mat + sigma_e_val3^2 + sigma_epsilon_val3^2) ) 
  #wrong.....
  #I think I need 2 steps...
  #2 emulators: 
  #     including the points to find expected values
  #     excluding the points to find variances
  #step 1) find expected values of proposed points
  #step 2) calculate variance INCLUDING the new points in emulator
  #step 3) calculate implausibility using 
  #     (i) expected values EXCLUDING new points
  #     (ii) variances INCLUDING new points
  #end of wrong.
  
  
  
  #redo above:
  
  #Aim: v90 of non-implausible region
  #wait wait wait 
  #non-imp stays same!!
  #but new variances (over same region) calculated using new emulator including the points
  
  #1)steps previously done to find non-implausible region
  #2)but then use these points as xP in NEW emulator 
  #using D=f(xD) because finding variance (doesn't depend on D)
  #then extract v90 (which by setup will be at each xP point)
  #then calculate v90 of extraction using quantile!
  
  #part 1
  xD_1st <- matrix(oldpts_vec_val3,ncol=2,byrow = TRUE)
  D_1st <- val3_actual_stellarmassfunctions
  x_grid <- seq(-0.001,1.001,len=grid_length_val3)
  xP_1st <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP_1st,xD=xD_1st,D=D_1st,theta=theta_val3,sigma=sigma_val3,E_f=E_f_val3,nugget=nugget_val3)
  E_D_fx <- em_out[,"ExpD_f(x)"]
  Var_D_fx <- em_out[,"VarD_f(x)"]
  
  Imp <- sqrt( (E_D_fx - z_val3)^2 / (Var_D_fx + sigma_e_val3^2 + sigma_epsilon_val3^2) ) 
  
  nonimplausible_indices <- which(Imp<3)
  nonimplausible_points <- xP_1st[nonimplausible_indices,]
  
  
  #part 2
  xD_2nd <- matrix(c(oldpts_vec_val3,newpts_vec_val3),ncol=2,byrow=TRUE)
  D_2nd <- f(xD_2nd)
  xP_2nd <- nonimplausible_points
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP_2nd,xD=xD_2nd,D=D_2nd,theta=theta_val3,sigma=sigma_val3,E_f=E_f_val3,nugget=nugget_val3)
  Var_D <- em_out[,"VarD_f(x)"]
  
  v90_W2 <- quantile(Var_D,percentile_val3)
  return(v90_W2)
}

M_function <- function(newpts_val, 
                       grid_length_val3,
                       xD_val3,
                       theta_val3,
                       sigma_val3,
                       E_f_val3,
                       nugget_val3,
                       
                       z_val3,
                       sigma_e_val3,
                       sigma_epsilon_val3,
                       
                       val3_actual_stellarmassfunctions){
  
  xD <- xD_val3
  D <- val3_actual_stellarmassfunctions
  
  #this is the key change: xP only the new proposed points!!!!!
  xP <- matrix(newpts_val,ncol=2,byrow=TRUE)
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=theta_val3,sigma=sigma_val3,E_f=E_f_val3,nugget=nugget_val3)  
  
  #change this slightly too
  E_D_fx <- em_out[,"ExpD_f(x)"] 
  Var_D_fx <- em_out[,"VarD_f(x)"]
  
  Imp <- sqrt( (E_D_fx - z_val3)^2 / (Var_D_fx + sigma_e_val3^2 + sigma_epsilon_val3^2) ) 
  
  M <- max(Imp)
  
  return(M)
}

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

M_sigmoidal <- function(M_val){
  #we choose b and d the same as Tc sigmoidal
  #a is steeper so approaches at 2.9 reaching middle at 3
  #b is centred at 3 (implausibility boundary)
  a <- 50 #steepness
  b <- 3 #height
  c <- 3 #where centred (x value halfway up the sigmoidal)
  d <- 0.5 #steepness of linear term added on
  x <- M_val
  if (M_val<=c){
    sigmoidal <- b/(1+exp(-a*(x-c)))
  }
  if (M_val>c){
    sigmoidal <- b/(1+exp(-a*(x-c))) + d*x - d*c
  }
  return(sigmoidal)
}


#this is new version of "sigmoidal_optim_function"
big_utility_function <- function(newpts, #will optimise first argument
                                 oldpts,
                                 grid_length=50,
                                 theta_val2=0.5,
                                 sigma_val2=0.5,
                                 E_f_val2=0,
                                 nugget_val2=10^(-6),
                                 true_runs_vec,
                                 
                                 ceiling,
                                 
                                 runtime_emulator_theta,
                                 runtime_emulator_sigma,
                                 runtime_emulator_E,
                                 runtime_emulator_nugget,
                                 
                                 percentile_val2,
                                 
                                 #new for this function:
                                 z_val2,
                                 sigma_e_val2,
                                 sigma_epsilon_val2,
                                 
                                 val2_actual_runtimes,
                                 val2_actual_stellarmassfunctions
){
  
  oldpts_matrix <- matrix(oldpts,ncol=2,byrow=TRUE) 
  newpts_matrix <- matrix(newpts,ncol=2,byrow=TRUE)
  xD <- rbind(oldpts_matrix,newpts_matrix)
  
  
  #V_{90}^{W2}
  v90_W2 <- V90W2_function(grid_length_val3=grid_length,
                           
                           oldpts_vec_val3=oldpts,
                           newpts_vec_val3=newpts,
                           
                           theta_val3=theta_val2,
                           sigma_val3=sigma_val2,
                           E_f_val3=E_f_val2,
                           nugget_val3=nugget_val2,
                           
                           z_val3=z_val2,
                           sigma_e_val3=sigma_e_val2,
                           sigma_epsilon_val3=sigma_epsilon_val2,
                           
                           percentile_val3=percentile_val2,
                           
                           val3_actual_stellarmassfunctions=val2_actual_stellarmassfunctions)
  
  #sigmoid(Tc)
  Tc <- emulate_log_run_time_Tc(new_batch_vector=newpts,
                                xD_matrix=oldpts_matrix,
                                D_vector=val2_actual_runtimes,
                                theta_val=runtime_emulator_theta,
                                sigma_val=runtime_emulator_sigma,
                                E_f_val=runtime_emulator_E,
                                nugget_val=runtime_emulator_nugget)
  
  sigmoid_Tc <- Tc_sigmoidal(Tc_val = Tc,
                             ceiling_val = ceiling)
  
  #sigmoid(M)
  M <- M_function(newpts_val=newpts, 
                  grid_length_val3=grid_length,
                  xD_val3=oldpts_matrix,
                  theta_val3=theta_val2,
                  sigma_val3=sigma_val2,
                  E_f_val3=E_f_val2,
                  nugget_val3=nugget_val2,
                  
                  z_val3=z_val2,
                  sigma_e_val3=sigma_e_val2,
                  sigma_epsilon_val3=sigma_epsilon_val2,
                  
                  val3_actual_stellarmassfunctions=val2_actual_stellarmassfunctions)
  
  sigmoid_M <- M_sigmoidal(M)
  
  
  
  final <- 1/2*log10(v90_W2) + sigmoid_Tc + sigmoid_M
  return(final)
}


#to create graph of history matching
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
                    xD=rbind(xD,matrix(hm_new_vec,ncol=2,byrow = TRUE)),
                    x_grid=x_grid,
                    xD_col=rep(c("purple","pink"),c(nrow(xD),nrow(matrix(hm_new_vec,ncol=2,byrow = TRUE)))),
                    color.palette=imp_cols,main="Implausibility I(x) and Best Design")
  
}


#I will try and plot the variance in just the non-implausible region now...
hm_var_plot <- function(xD_vec_original,
                        xD_vec_new,
                        grid_length, 
                        theta_val, 
                        sigma_val, 
                        E_f_val,
                        
                        nugget_val,
                        
                        hm_var_z,
                        hm_var_sigma_e,
                        hm_var_sigma_epsilon,
                        
                        hm_D){
  
  #again!!! need to do it in two steps...
  #1) calc implausibility
  #2) calc new var and put na's where implausible
  
  par(mar=c(4,4,3,3))
  var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
  
  #1)
  xD_1st <- matrix(xD_vec_original,ncol=2,byrow = TRUE)
  D_1st <- hm_D
  
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD_1st,D=D_1st,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)   
  
  E_D_fx_mat <- em_out[,"ExpD_f(x)"] 
  Var_D_fx_mat <- em_out[,"VarD_f(x)"]
  
  Imp_mat <- sqrt( (E_D_fx_mat - hm_var_z)^2 / (Var_D_fx_mat + hm_var_sigma_e^2 + hm_var_sigma_epsilon^2) ) 
  
  nonimplausible_indices <- which(Imp_mat<3)
  
  
  
  #2)
  xD_2nd <- matrix(c(xD_vec_original,xD_vec_new),ncol=2,byrow=TRUE)
  D_2nd <- f(xD_2nd)
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD_2nd,D=D_2nd,theta=theta_val,sigma=sigma_val,E_f=E_f_val,nugget=nugget_val)   
  
  all_var <- em_out[,"VarD_f(x)"]
  all_var[-nonimplausible_indices] <- NA
  
  Var_D_fx_mat <- matrix(all_var,nrow=length(x_grid),ncol=length(x_grid)) 
  
  emul_fill_cont(cont_mat=Var_D_fx_mat,cont_levs=NULL,xD=xD_2nd,x_grid=x_grid,
                 color.palette=var_cols,main="Emulator Adjusted Variance Var_D[f(x)]",
                 #new part:
                 xD_col=rep(c("purple","pink"),c(length(xD_vec_original)/2,length(xD_vec_new)/2)) #/2 because x and y coords in vector so number of points is length/2
  )
  
  
}



#history matching changed big function
hm_everything_function_sigmoidal <- function(want_all_var_graphs,
                                             want_only_best_var_graph,
                                             want_freeze, #initialise design at melt design found now with harsher sigmoidal
                                             
                                             batch0,
                                             #batch0_evaluated, #evaluated run times
                                             
                                             actual_runtimes,
                                             actual_stellarmassfunctions,
                                             
                                             
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
                                             
                                             ceiling_value, #NOT LOGGED max total run time cannot exceed
                                             
                                             
                                             #NEW FOR HM:
                                             z_observation,
                                             sigma_e,
                                             sigma_epsilon
){
  
  #create random function (won't affect variance but function needs it)
  f <- function(x) -sin(2*pi*x[,2]) + 0.9*sin(2*pi*(1-x[,1])*(1-x[,2]))
  
  
  #must redo this and place initial points in the non-implausible region
  
  # #determine initial starting points 
  # #dependent on batch_size (number of new points) we want
  # if (batch_size==4){
  #   starting_pts <- c(0.9,0.85,
  #                     0.6,0.95, 
  #                     0.4,0.05,
  #                     0.1,0.35)
  # }
  # if (batch_size==6){
  #   starting_pts <- c(0.6,0.95,
  #                     0.9,0.85,
  #                     0.95,0.65,
  #                     
  #                     0.05,0.4,
  #                     0.1,0.15,
  #                     0.35,0.05)
  # }
  # if (batch_size==8){
  #   starting_pts <- c(0.6,0.95,
  #                     0.9,0.85,
  #                     0.95,0.65,
  #                     
  #                     0.05,0.4,
  #                     0.1,0.15,
  #                     0.35,0.05,
  #                     
  #                     0.3,0.5,
  #                     0.7,0.6)
  # }
  # if (batch_size==10){
  #   starting_pts <- c(0.6,0.95,
  #                     0.9,0.85,
  #                     0.95,0.65,
  #                     
  #                     0.05,0.4,
  #                     0.1,0.15,
  #                     0.35,0.05,
  #                     
  #                     0.3,0.5,
  #                     0.7,0.6,
  #                     
  #                     0.2,0.8,
  #                     0.8,0.1)
  # }
  # if (batch_size==5){
  #   starting_pts <- c(0.6,0.95,
  #                     0.9,0.85,
  #                     0.95,0.65,
  #                     
  #                     0.4,0.05,
  #                     0.1,0.35)
  # }
  # if (batch_size==7){
  #   starting_pts <- c(0.6,0.95,
  #                     0.9,0.85,
  #                     0.95,0.65,
  #                     
  #                     0.05,0.4,
  #                     0.1,0.15,
  #                     0.35,0.05,
  #                     
  #                     0.7,0.6)
  #   
  # }
  # if (batch_size==9){
  #   starting_pts <- c(0.6,0.95,
  #                     0.9,0.85,
  #                     0.95,0.65,
  #                     
  #                     0.05,0.4,
  #                     0.1,0.15,
  #                     0.35,0.05,
  #                     
  #                     0.3,0.5,
  #                     0.7,0.6,
  #                     
  #                     0.2,0.8)
  # }
  # 
  #this creates "number_ICs" lots of initial starting points
  #1st is with no noise
  #2nd onwards are with 2D multivariate normal noise mean (0,0) var=(sigma_noise^2 , 0
  #                                                                  0 , sigma_noise^2)
  # noise_vector <- noise_function(batch_size,
  #                                number_ICs,
  #                                sigma_noise)
  # 
  # starting_pts_without_ICs <- rep(starting_pts,number_ICs)
  # starting_pts_with_ICs <- starting_pts_without_ICs+noise_vector
  # 
  # 
  # 
  # 
  # #this ensures any noisy initial values below 0 are changed to 0.01
  # #similarly above 1 changed to 0.99
  # starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs<0,0.01)
  # starting_pts_with_ICs <- replace(starting_pts_with_ICs,starting_pts_with_ICs>1,0.99)
  # 
  
  xD <- matrix(batch0,ncol=2,byrow = TRUE)
  D <- f(xD)
  
  x_grid <- seq(-0.001,1.001,len=grid_length)
  xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))
  
  em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=actual_stellarmassfunctions,theta=theta_emulator,sigma=sigma_emulator,E_f=E_emulator,nugget=nugget_emulator)  
  E_D_fx_mat <- em_out[,"ExpD_f(x)"] 
  Var_D_fx_mat <- em_out[,"VarD_f(x)"]
  
  Imp_mat <- sqrt( (E_D_fx_mat - z_observation)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2) ) 
  
  nonimplausible_indices <- which(Imp_mat<3)
  non_implausible_points <- xP[nonimplausible_indices,]
  
  list_starting <- list()
  number_nonimp <- nrow(non_implausible_points)
  for (i in 1:number_ICs){
    random_indices <- sample(1:number_nonimp,size=batch_size)
    list_starting[[i]] <- c(t(non_implausible_points[random_indices,]))
  }
  
  starting_pts_with_ICs <- unlist(list_starting)
  
  
  
  #need to change this
  
  # #now we find the optim value over all ICs
  # #we save the parameters and the optimised value
  # optim_par <- rep(0,number_ICs*batch_size*2)
  # optim_value <- rep(0,number_ICs)
  # for (j in 1:number_ICs){
  #   optim <- optim(method = "L-BFGS-B",
  #                  
  #                  par=starting_pts_with_ICs[(2*(j-1)*batch_size+1):(2*j*batch_size)], 
  #                  
  #                  fn=sigmoidal_optim_function,
  #                  
  #                  oldpts=batch0,
  #                  grid_length=grid_length,
  #                  theta_val2=theta_emulator,
  #                  sigma_val2=sigma_emulator,
  #                  E_f_val2=E_emulator,
  #                  nugget_val2=nugget_emulator,
  #                  true_runs_vec=batch0_evaluated,
  #                  
  #                  runtime_emulator_theta=run_time_emulator_theta,
  #                  runtime_emulator_sigma=run_time_emulator_sigma,
  #                  runtime_emulator_E=run_time_emulator_E,
  #                  runtime_emulator_nugget=run_time_emulator_nugget,
  #                  
  #                  ceiling=ceiling_value,
  #                  freeze=FALSE, #this step is just about doing the melt (gentle) sigmoidal
  #                  
  #                  lower=0,
  #                  upper=1)
  #   optim_par[((batch_size*2)*(j-1)+1):((batch_size*2)*j)] <- optim$par
  #   optim_value[j] <- optim$value
  # }
  # 
  
  #now we find the optim value over all ICs
  #we save the parameters and the optimised value
  optim_par <- rep(0,number_ICs*batch_size*2)
  optim_value <- rep(0,number_ICs)
  for (j in 1:number_ICs){
    optim <- optim(method = "L-BFGS-B",
                   
                   par=starting_pts_with_ICs[(2*(j-1)*batch_size+1):(2*j*batch_size)],
                   
                   fn=big_utility_function,
                   
                   oldpts=batch0,
                   grid_length=grid_length,
                   theta_val2=theta_emulator,
                   sigma_val2=sigma_emulator,
                   E_f_val2=E_emulator,
                   nugget_val2=nugget_emulator,
                   
                   #true_runs_vec=batch0_evaluated,
                   val2_actual_runtimes=actual_runtimes,
                   val2_actual_stellarmassfunctions=actual_stellarmassfunctions,
                   
                   
                   runtime_emulator_theta=run_time_emulator_theta,
                   runtime_emulator_sigma=run_time_emulator_sigma,
                   runtime_emulator_E=run_time_emulator_E,
                   runtime_emulator_nugget=run_time_emulator_nugget,
                   
                   ceiling=ceiling_value,
                   
                   percentile_val2 = percentile,
                   
                   #NEW:
                   z_val2=z_observation,
                   sigma_e_val2=sigma_e,
                   sigma_epsilon_val2=sigma_epsilon,
                   
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
  
  
  #need to check over this...
  #return M as well...
  
  #we calculate v90 and Tc for each set of optimised locations
  #calculate V90 and Tc for each set of locations
  oldpts_matrix <- matrix(batch0,ncol=2,byrow=TRUE) 
  v90_W2 <- rep(0,number_ICs)
  Tc <- rep(0,number_ICs)
  M <- rep(0,number_ICs)
  for (z in 1:(number_ICs)){ 
    newpts_matrix <- matrix(optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],ncol=2,byrow=TRUE)
    xD <- rbind(oldpts_matrix,newpts_matrix)
    D <- f(xD)
    
    v90_W2[z] <- V90W2_function(grid_length_val3=grid_length,
                                
                                oldpts_vec_val3=c(t(oldpts_matrix)),
                                newpts_vec_val3=c(t(newpts_matrix)),
                                
                                theta_val3=theta_emulator,
                                sigma_val3=sigma_emulator,
                                E_f_val3=E_emulator,
                                nugget_val3=nugget_emulator,
                                
                                z_val=z_observation,
                                sigma_e_val=sigma_e,
                                sigma_epsilon_val=sigma_epsilon,
                                
                                percentile_val3=percentile,
                                
                                val3_actual_stellarmassfunctions=actual_stellarmassfunctions)
    
    Tc[z] <- emulate_log_run_time_Tc(new_batch_vector=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)],
                                     xD_matrix=oldpts_matrix,
                                     D_vector=actual_runtimes,
                                     theta_val=run_time_emulator_theta,
                                     sigma_val=run_time_emulator_sigma,
                                     E_f_val=run_time_emulator_E,
                                     nugget_val=run_time_emulator_nugget)
    
    M[z] <- M_function(newpts_val=optim_par[((batch_size*2)*(z-1)+1):((batch_size*2)*z)], 
                       grid_length_val3=grid_length,
                       xD_val3=oldpts_matrix,
                       theta_val3=theta_emulator,
                       sigma_val3=sigma_emulator,
                       E_f_val3=E_emulator,
                       nugget_val3=nugget_emulator,
                       
                       z_val3=z_observation,
                       sigma_e_val3=sigma_e,
                       sigma_epsilon_val3=sigma_epsilon,
                       
                       val3_actual_stellarmassfunctions=actual_stellarmassfunctions)
  }
  
  
  
  
  
  #best design (over the ICs) based on which gives the lowest optim value
  minimum_index <- which.min(optim_value)
  
  best_design_v90_W2 <- v90_W2[minimum_index]
  best_design_Tc <- Tc[minimum_index]
  best_design_M <- M[minimum_index]
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
  
  hm_plot(hm_z=z_observation,
          hm_sigma_e=sigma_e,
          hm_sigma_epsilon=sigma_epsilon,
          hm_grid_length=grid_length,
          hm_old_vec=batch0,
          hm_new_vec=c(t(best_design_locations)),
          hm_D_runs=actual_stellarmassfunctions,
          hm_theta=theta_emulator,
          hm_sigma=sigma_emulator,
          hm_E_f=E_emulator,
          hm_nugget=nugget_emulator)
  
  hm_var_plot(xD_vec_original=batch0,
              xD_vec_new=c(t(best_design_locations)),
              grid_length=grid_length, 
              theta_val=theta_emulator, 
              sigma_val=sigma_emulator, 
              E_f_val=E_emulator,
              nugget_val=nugget_emulator,
              
              hm_var_z=z_observation,
              hm_var_sigma_e=sigma_e,
              hm_var_sigma_epsilon=sigma_epsilon,
              
              hm_D=actual_stellarmassfunctions)
  
  #return some things we may want to look at
  return(list("All v90 W2"=v90_W2,
              "All Tc"=Tc,
              "All M"=M,
              "All locations"=optim_par,
              "All values"=optim_value,
              
              "Best Design v90 W2"=best_design_v90_W2,
              "Best Design Tc"=best_design_Tc,
              "Best Design M"=best_design_M,
              "Best Design locations"=best_design_locations,
              "Best Design value"=best_design_value,
              
              "Freeze v90"=freeze_v90,
              "Freeze Tc"=freeze_Tc,
              "Freeze location"=freeze_locations,
              "Freeze value"=freeze_value,
              
              "Best design index"=minimum_index
  ))
}



## now run this  -----------------------------------------------------------


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





## this is code needed for xD ----------------------------------------------
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


## this is code needed for D runtimes --------------------------------------

D_runtimes <- log10(matrix_sim_run_times[,"Runtime [CPU-Hours]"])
sd_runtimes <- sd(D_runtimes)
mean_runtimes <- mean(D_runtimes)
nugget_runtimes <- 0.145064^2

## 9: this is code needed for D, hm ----------------------------------------
D_9 <- c(0.079975008,
         0.021868166,
         0.002343018,
         0.070759138,
         0.044673540,
         0.005935645,
         0.059668853,
         0.054201812,
         0.012652296)
sd_9 <- sd(D_9)
mean_9 <- mean(D_9)
nugget_9 <- 0.05430212^2
z_9 <- 0.01225
#CHANGING ERRORS !!!!
sigma_e_9 <- 0.001225
sigma_epsilon_9 <- 0.001225

## 10.5: this is code needed for D, hm -------------------------------------
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
z_10.5 <- 0.00512
#CHANGING ERRORS !!!!
sigma_e_10.5 <- 0.000512
sigma_epsilon_10.5 <- 0.000512




## run this to find the optimal 6-point design large ceiling ---------------


# png(file="pts6_10000_%d.png",width=1400,height=1200,res=200)
# hm_pts6_10000 <- hm_everything_function_sigmoidal(want_all_var_graphs=FALSE,
#                                                   want_only_best_var_graph=TRUE,
#                                                   want_freeze=FALSE, #initialise design at melt design found now with harsher sigmoidal
# 
#                                                   batch0=c(t(xD)),
#                                                   #batch0_evaluated=D, #evaluated run times
# 
#                                                   actual_runtimes=D_runtimes,
#                                                   actual_stellarmassfunctions=D_10.5,
# 
#                                                   grid_length=60,
# 
#                                                   batch_size=6, #4, 6, 8 or 10
# 
#                                                   theta_emulator=1/2,
#                                                   sigma_emulator=sd_10.5,
#                                                   E_emulator=mean_10.5,
#                                                   nugget_emulator=nugget_10.5, #this is the nugget of hm 9!!!
# 
#                                                   run_time_emulator_theta=1,
#                                                   run_time_emulator_sigma=sd_runtimes,
#                                                   run_time_emulator_E=mean_runtimes,
#                                                   run_time_emulator_nugget=nugget_runtimes,
# 
#                                                   percentile=0.9, #percentile of variance used in optim
#                                                   number_ICs=40, #how many different starting points we want
#                                                   sigma_noise=0.15, #N(0,sigma_noise^2) for random noise
# 
#                                                   ceiling_value=10000, #NOT LOGGED max total run time cannot exceed
# 
# 
#                                                   #NEW:
#                                                   z_observation=z_10.5,
#                                                   sigma_e=sigma_e_10.5,
#                                                   sigma_epsilon=sigma_epsilon_10.5
# )
# dev.off()
# 
# hm_pts6_10000$`Best Design v90 W2`
# #4.333608e-07
# hm_pts6_10000$`Best Design Tc`
# #5362.355
# hm_pts6_10000$`Best Design M`
# #2.563278
# hm_pts6_10000$`Best Design locations`
# # [,1]       [,2]
# # [1,] 0.9999822 0.34746227
# # [2,] 0.9083349 0.95194057
# # [3,] 0.9967341 0.73509330
# # [4,] 0.2457146 0.03650289
# # [5,] 0.4968689 0.06203308
# # [6,] 0.5361542 0.99974680
# hm_pts6_10000$`Best Design value`
# #-3.181575





# Figure 30 ---------------------------------------------------------------



#initial idea, may change as I go...
#1) Calculate log10 run time emulator
#2) Extract emulator expectation and emulator variance at each point of new design
#         i.e. xP = new design
#3) At each point:
#         (i) simulate 1000 runs from Normal distribution with mean emulator expectation and variance emulator variance
#         (ii) 10^ of each realisation
#4) tc_total <- rep(0,1000)
#5) For i in 1 to 1000:
#         (i) tc_total[i] = sum of ith simulation at each point
#6) Calculate mean(tc_total)
#...
#7) Plot histogram of Tc,totals
#8) Put a vertical line of Tc,total using emulator expectations only
#9) Put a vertical line of mean(tc_total)
#10) Put vertical lines at quantile(tc_total,c(0.025,0.975))


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

locations <- c(0.9999822, 0.34746227,
               0.9083349, 0.95194057,
               0.9967341, 0.73509330,
               0.2457146, 0.03650289,
               0.4968689, 0.06203308,
               0.5361542, 0.99974680)

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

D_runtimes <- log10(matrix_sim_run_times[,"Runtime [CPU-Hours]"])
sd_runtimes <- sd(D_runtimes)
mean_runtimes <- mean(D_runtimes)
nugget_runtimes <- 0.145064^2

library(pdist)

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




results <- simulation_Tc(sim_new_pts_vec=locations,
                         sim_old_pts_vec=c(t(xD)),
                         sim_runtime_D=D_runtimes,
                         sim_runtime_theta=1,
                         sim_runtime_sigma=sd_runtimes,
                         sim_runtime_E=mean_runtimes,
                         sim_runtime_nugget=nugget_runtimes,
                         sim_batchsize=6)


#png(file="simulated_histogram_final.png",width=2000,height=1200,res=200)
hist(results$`simulated totals`,
     xlab="Total Run Time (natural scale)",
     main = "Histogram of 1000 Simulated Total Run Times")
abline(v=c(results$`expected total`,results$`mean total`,quantile(results$`simulated totals`,c(0.025,0.975))),
       col=c("red","gold","green","blue"),
       lty=c(1,2,1,1),
       lwd=2.5)
legend(x="topright",
       legend=c("Expectation","Mean","2.5th Percentile","97.5th Percentile"),
       lty=c(1,2,1,1),
       col=c("red","gold","green","blue"),
       lwd=2.5)
#dev.off()

quantile(results$`simulated totals`,c(0.025,0.975))
# 2.5%      97.5% 
# 4935.638  5876.317 

results$`expected total` 
#5362.355

results$`mean total`
#5380.732

# saveRDS(results, file="results.RData")

# saveRDS(results, file="results2.RData")





# Figure 31 ---------------------------------------------------------------



## 31a ---------------------------------------------------------------------

# See the 6-point design with 6000 ceiling
# Figure 23c


## 31b ---------------------------------------------------------------------



### these are the new functions ---------------------------------------------


#key idea:
#We want E[U] = E[ -(  1/2*log10(v90) + sigmoid(Tc) ) ]
# = -(  1/2*log10(v90) + E[ sigmoid(Tc) ] )
# ≈ -(  1/2*log10(v90) + 1/1000 \sum{i=1}{1000} sigmoid(Tc^(i))  )
#in other words...
#calculate v90 as normal
#but instead of using Tc
#calculate all the simulated Tcs
#plug each one into sigmoidal function
#then put 1/1000 sum of those sigmoidal Tc's into the criteria!


#now I will incorporate this criteria into the original criteria I used (pre hm)
#then I will use the new criteria to find a Wave 1 design using 6 points
#see if more risk averse

#this is fn I created in task 2:
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


### this contains all functions needed (including new ones above) --------------------------------------------------------------------



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





# design ------------------------------------------------------------------


# png(file="simulations_pts6_6000_%d.png",width=1400,height=1200,res=200)
# simulations_pts6_6000 <- new_simulation_everything_function_sigmoidal(want_all_var_graphs=FALSE,
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
#                                                                       ceiling_value=6000 #NOT LOGGED max total run time cannot exceed
# )
# dev.off()
# saveRDS(simulations_pts6_6000, file="simulations_pts6_6000.RData")









# End of Section 9 --------------------------------------------------------
