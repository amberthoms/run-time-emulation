
# Section 8: History Matching ---------------------------------------------

# Figures 25-28

# Functions ---------------------------------------------------------------

# Figure 25 ---------------------------------------------------------------


observed_data_matrix <- matrix(ncol = 5,
                               byrow = TRUE,
                               data = c(
                                 6.25, 0.50, 31.1, 21.6, 9,
                                 6.75, 0.50, 18.1, 6.6, 19,
                                 7.10, 0.20, 17.9, 5.7, 18,
                                 7.30, 0.20, 43.1, 8.7, 46,
                                 7.50, 0.20, 31.6, 9.0, 51,
                                 7.70, 0.20, 34.8, 8.4, 88,
                                 7.90, 0.20, 27.3, 4.2, 140,
                                 8.10, 0.20, 28.3, 2.8, 243,
                                 8.30, 0.20, 23.5, 3.0, 282,
                                 8.50, 0.20, 19.2, 1.2, 399,
                                 8.70, 0.20, 18.0, 2.6, 494,
                                 8.90, 0.20, 14.3, 1.7, 505,
                                 9.10, 0.20, 10.2, 0.6, 449,
                                 9.30, 0.20, 9.59, 0.55, 423,
                                 9.50, 0.20, 7.42, 0.41, 340,
                                 9.70, 0.20, 6.21, 0.37, 290,
                                 9.90, 0.20, 5.71, 0.35, 268,
                                 10.10, 0.20, 5.51, 0.34, 260,
                                 10.30, 0.20, 5.48, 0.34, 259,
                                 10.50, 0.20, 5.12, 0.33, 242,
                                 10.70, 0.20, 3.55, 0.27, 168,
                                 10.90, 0.20, 2.41, 0.23, 114,
                                 11.10, 0.20, 1.27, 0.16, 60,
                                 11.30, 0.20, 0.338, 0.085, 16,
                                 11.50, 0.20, 0.042, 0.030, 2,
                                 11.70, 0.20, 0.021, 0.021, 1,
                                 11.90, 0.20, 0.042, 0.030, 2
                               )
)

colnames(observed_data_matrix) <- c("mid point",
                                    "bin width",
                                    "phi",
                                    "error",
                                    "number")

View(observed_data_matrix)





Run_0 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/Run_0/halo_catalogue.txt")
Run_1 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/Run_1/halo_catalogue.txt")
Run_2 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/Run_2/halo_catalogue.txt")
Run_3 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/Run_3/halo_catalogue.txt")
Run_4 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/Run_4/halo_catalogue.txt")
Run_5 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/Run_5/halo_catalogue.txt")
Run_6 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/Run_6/halo_catalogue.txt")
Run_7 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/Run_7/halo_catalogue.txt")
Run_8 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/Run_8/halo_catalogue.txt")

#we add 10 to log10(...) which is equivalent to log10(...x10^10)
hist_run_0 <- hist(log10(Run_0$V2)+10, 
                   plot=F,
                   breaks=seq(6.5,12,by=0.5))

hist_run_1 <- hist(log10(Run_1$V2)+10, 
                   plot=F,
                   breaks=seq(6.5,12,by=0.5))

hist_run_2 <- hist(log10(Run_2$V2)+10, 
                   plot=F,
                   breaks=seq(6.5,12,by=0.5))

hist_run_3 <- hist(log10(Run_3$V2)+10, 
                   plot=F,
                   breaks=seq(6.5,12,by=0.5))

hist_run_4 <- hist(log10(Run_4$V2)+10, 
                   plot=F,
                   breaks=seq(6.5,12,by=0.5))

hist_run_5 <- hist(log10(Run_5$V2)+10, 
                   plot=F,
                   breaks=seq(6.5,12,by=0.5))

hist_run_6 <- hist(log10(Run_6$V2)+10, 
                   plot=F,
                   breaks=seq(6.5,12,by=0.5))

hist_run_7 <- hist(log10(Run_7$V2)+10, 
                   plot=F,
                   breaks=seq(6.5,12,by=0.5))

hist_run_8 <- hist(log10(Run_8$V2)+10, 
                   plot=F,
                   breaks=seq(6.5,12,by=0.5))


#number of galaxies per log_10 (M_st)
run0_step1 <- hist_run_0$counts/0.5
run1_step1 <- hist_run_1$counts/0.5
run2_step1 <- hist_run_2$counts/0.5
run3_step1 <- hist_run_3$counts/0.5
run4_step1 <- hist_run_4$counts/0.5
run5_step1 <- hist_run_5$counts/0.5
run6_step1 <- hist_run_6$counts/0.5
run7_step1 <- hist_run_7$counts/0.5
run8_step1 <- hist_run_8$counts/0.5


#number of galaxies per log10(M_st) per unit volume
run0_step2 <- run0_step1/6402
run1_step2 <- run1_step1/6402
run2_step2 <- run2_step1/6402
run3_step2 <- run3_step1/6402
run4_step2 <- run4_step1/6402
run5_step2 <- run5_step1/6402
run6_step2 <- run6_step1/6402
run7_step2 <- run7_step1/6402
run8_step2 <- run8_step1/6402



#volume averaged stellar mass function
library(ggplot2)
library(viridis)
data <- data.frame(
  x_for_all=hist_run_0$mids,
  run0=run0_step2,
  run1=run1_step2,
  run2=run2_step2,
  run3=run3_step2,
  run4=run4_step2,
  run5=run5_step2,
  run6=run6_step2,
  run7=run7_step2,
  run8=run8_step2)

#need to multiply table phi by 10^(-3)
observed <- data.frame(
  observed_data_matrix_x = observed_data_matrix[,"mid point"],
  observed_data_matrix_y = 10^(-3)*observed_data_matrix[,"phi"]
)



observed <- data.frame(
  observed_data_matrix_x = observed_data_matrix[,"mid point"],
  observed_data_matrix_y = 10^(-3)*observed_data_matrix[,"phi"],
  errors = 10^(-3)*observed_data_matrix[,"error"]
)

myerrors <- data.frame(
  x_bins = c(9,10.5),
  my_errors = c(sqrt(0.001225^2+0.001225^2),sqrt(0.000512^2+0.000512^2)),
  y_vals = c(0.01225,0.00512)
)

#graph without stochastic repetitions
ggplot(data, aes(x=x_for_all)) +
  geom_line(aes(y=run0, colour = "run 0")) +
  geom_line(aes(y=run1, colour = "run 1")) +
  geom_line(aes(y=run2, colour = "run 2")) +
  geom_line(aes(y=run3, colour = "run 3")) +
  geom_line(aes(y=run4, colour = "run 4")) +
  geom_line(aes(y=run5, colour = "run 5")) +
  geom_line(aes(y=run6, colour = "run 6")) +
  geom_line(aes(y=run7, colour = "run 7")) +
  geom_line(aes(y=run8, colour = "run 8")) +
  geom_line(data=observed,aes(x=observed_data_matrix_x,y=observed_data_matrix_y, colour = "observed")) +
  ggtitle("Volume Averaged log Stellar Mass Functions") +
  xlab("log Galaxy Mass") +
  ylab("log Galaxy Counts") +
  scale_y_continuous(trans='log10') +
  #xlim(8.5,12) +
  geom_errorbar(data=observed, 
                aes(x=observed_data_matrix_x,ymin=pmax(observed_data_matrix_y - 2*errors, 0), ymax=observed_data_matrix_y+2*errors), 
                width=.05,
                position=position_dodge(0.05),
                color = "black") +
  geom_point(data=observed, 
             aes(x=observed_data_matrix_x,y=observed_data_matrix_y),
             color = "black")+
  geom_errorbar(data=myerrors, 
                aes(x=x_bins,ymin=y_vals-2*my_errors, ymax=y_vals+2*my_errors), 
                width=.05,
                position=position_dodge(0.05),
                color = "brown") +
  coord_cartesian(xlim =c(8.5, 12),ylim=c(0.00001,NA))+
  scale_colour_manual("", 
                      breaks = c("observed","run 0","run 1","run 2","run 3","run 4","run 5","run 6","run 7","run 8"),
                      values = c("observed"="black","run 0"="red","run 1"="orange","run 2"="gold","run 3"="green","run 4"="skyblue1","run 5"="royalblue3","run 6"="purple","run 7"="hotpink","run 8"="pink"))




#redone with units of axis

#png(file="supernew_stellar_plot.png",width=1760,height=1320,res=200)

ggplot(data, aes(x=x_for_all)) +
  geom_line(aes(y=run0, colour = "run 0")) +
  geom_line(aes(y=run1, colour = "run 1")) +
  geom_line(aes(y=run2, colour = "run 2")) +
  geom_line(aes(y=run3, colour = "run 3")) +
  geom_line(aes(y=run4, colour = "run 4")) +
  geom_line(aes(y=run5, colour = "run 5")) +
  geom_line(aes(y=run6, colour = "run 6")) +
  geom_line(aes(y=run7, colour = "run 7")) +
  geom_line(aes(y=run8, colour = "run 8")) +
  geom_line(data=observed,aes(x=observed_data_matrix_x,y=observed_data_matrix_y, colour = "observed")) +
  ggtitle("Volume Averaged log Stellar Mass Functions") +
  xlab(expression("log Galaxy Mass [M"[ʘ]*"]")) +
  ylab(expression("log Galaxy Counts ["*dex^{-1}~Mpc^{-3}*"]")) +
  scale_y_continuous(trans='log10') +
  #xlim(8.5,12) +
  geom_errorbar(data=observed, 
                aes(x=observed_data_matrix_x,ymin=pmax(observed_data_matrix_y - 2*errors, 0), ymax=observed_data_matrix_y+2*errors), 
                width=.05,
                position=position_dodge(0.05),
                color = "black") +
  geom_point(data=observed, 
             aes(x=observed_data_matrix_x,y=observed_data_matrix_y),
             color = "black")+
  geom_errorbar(data=myerrors, 
                aes(x=x_bins,ymin=y_vals-2*my_errors, ymax=y_vals+2*my_errors), 
                width=.05,
                position=position_dodge(0.05),
                color = "brown") +
  coord_cartesian(xlim =c(8.5, 12),ylim=c(0.00001,NA))+
  scale_colour_manual("", 
                      breaks = c("observed","run 0","run 1","run 2","run 3","run 4","run 5","run 6","run 7","run 8"),
                      values = c("observed"="black","run 0"="red","run 1"="orange","run 2"="gold","run 3"="green","run 4"="skyblue1","run 5"="royalblue3","run 6"="purple","run 7"="hotpink","run 8"="pink"))


#dev.off()








# Other: Investigating the stochastic repetitions and extracting outputs -------------------------
# Features graph with stochastic repetitions which is not featured in report
#1)make histogram for each
#2)make graph
#3)collect all 9 values
#4)collect all 10.5 values

stochastic_Run_0 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stochasticity_tests/Run_0/halo_catalogue.txt")
stochastic_Run_1 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stochasticity_tests/Run_1/halo_catalogue.txt")
stochastic_Run_2 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stochasticity_tests/Run_2/halo_catalogue.txt")
stochastic_Run_3 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stochasticity_tests/Run_3/halo_catalogue.txt")
stochastic_Run_4 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stochasticity_tests/Run_4/halo_catalogue.txt")
stochastic_Run_5 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stochasticity_tests/Run_5/halo_catalogue.txt")
stochastic_Run_6 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stochasticity_tests/Run_6/halo_catalogue.txt")
stochastic_Run_7 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stochasticity_tests/Run_7/halo_catalogue.txt")
stochastic_Run_8 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stochasticity_tests/Run_8/halo_catalogue.txt")
stochastic_Run_9 <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stochasticity_tests/Run_9/halo_catalogue.txt")

stochastic_hist_run_0 <- hist(log10(stochastic_Run_0$V2)+10, 
                              plot=F,
                              breaks=seq(6.5,12,by=0.5))

stochastic_hist_run_1 <- hist(log10(stochastic_Run_1$V2)+10, 
                              plot=F,
                              breaks=seq(6.5,12,by=0.5))

stochastic_hist_run_2 <- hist(log10(stochastic_Run_2$V2)+10, 
                              plot=F,
                              breaks=seq(6.5,12,by=0.5))

stochastic_hist_run_3 <- hist(log10(stochastic_Run_3$V2)+10, 
                              plot=F,
                              breaks=seq(6.5,12,by=0.5))

stochastic_hist_run_4 <- hist(log10(stochastic_Run_4$V2)+10, 
                              plot=F,
                              breaks=seq(6.5,12,by=0.5))

stochastic_hist_run_5 <- hist(log10(stochastic_Run_5$V2)+10, 
                              plot=F,
                              breaks=seq(6.5,12,by=0.5))

stochastic_hist_run_6 <- hist(log10(stochastic_Run_6$V2)+10, 
                              plot=F,
                              breaks=seq(6.5,12,by=0.5))

stochastic_hist_run_7 <- hist(log10(stochastic_Run_7$V2)+10, 
                              plot=F,
                              breaks=seq(6.5,12,by=0.5))

stochastic_hist_run_8 <- hist(log10(stochastic_Run_8$V2)+10, 
                              plot=F,
                              breaks=seq(6.5,12,by=0.5))

stochastic_hist_run_9 <- hist(log10(stochastic_Run_9$V2)+10, 
                              plot=F,
                              breaks=seq(6.5,12,by=0.5))




#number of galaxies per log_10 (M_st)
stochastic_run0_step1 <- stochastic_hist_run_0$counts/0.5
stochastic_run1_step1 <- stochastic_hist_run_1$counts/0.5
stochastic_run2_step1 <- stochastic_hist_run_2$counts/0.5
stochastic_run3_step1 <- stochastic_hist_run_3$counts/0.5
stochastic_run4_step1 <- stochastic_hist_run_4$counts/0.5
stochastic_run5_step1 <- stochastic_hist_run_5$counts/0.5
stochastic_run6_step1 <- stochastic_hist_run_6$counts/0.5
stochastic_run7_step1 <- stochastic_hist_run_7$counts/0.5
stochastic_run8_step1 <- stochastic_hist_run_8$counts/0.5
stochastic_run9_step1 <- stochastic_hist_run_9$counts/0.5


#number of galaxies per log10(M_st) per unit volume
stochastic_run0_step2 <- stochastic_run0_step1/6402
stochastic_run1_step2 <- stochastic_run1_step1/6402
stochastic_run2_step2 <- stochastic_run2_step1/6402
stochastic_run3_step2 <- stochastic_run3_step1/6402
stochastic_run4_step2 <- stochastic_run4_step1/6402
stochastic_run5_step2 <- stochastic_run5_step1/6402
stochastic_run6_step2 <- stochastic_run6_step1/6402
stochastic_run7_step2 <- stochastic_run7_step1/6402
stochastic_run8_step2 <- stochastic_run8_step1/6402
stochastic_run9_step2 <- stochastic_run9_step1/6402



#volume averaged stellar mass function
library(ggplot2)
library(viridis)
data <- data.frame(
  x_for_all=hist_run_0$mids,
  run0=run0_step2,
  run1=run1_step2,
  run2=run2_step2,
  run3=run3_step2,
  run4=run4_step2,
  run5=run5_step2,
  run6=run6_step2,
  run7=run7_step2,
  run8=run8_step2,
  stochastic_run0=stochastic_run0_step2,
  stochastic_run1=stochastic_run1_step2,
  stochastic_run2=stochastic_run2_step2,
  stochastic_run3=stochastic_run3_step2,
  stochastic_run4=stochastic_run4_step2,
  stochastic_run5=stochastic_run5_step2,
  stochastic_run6=stochastic_run6_step2,
  stochastic_run7=stochastic_run7_step2,
  stochastic_run8=stochastic_run8_step2,
  stochastic_run9=stochastic_run9_step2
)

#need to multiply table phi by 10^(-3)
observed <- data.frame(
  observed_data_matrix_x = observed_data_matrix[,"mid point"],
  observed_data_matrix_y = 10^(-3)*observed_data_matrix[,"phi"]
)

ggplot(data, aes(x=x_for_all)) +
  geom_line(aes(y=run0, colour = "run 0")) +
  geom_line(aes(y=run1, colour = "run 1")) +
  geom_line(aes(y=run2, colour = "run 2")) +
  geom_line(aes(y=run3, colour = "run 3")) +
  geom_line(aes(y=run4, colour = "run 4")) +
  geom_line(aes(y=run5, colour = "run 5")) +
  geom_line(aes(y=run6, colour = "run 6")) +
  geom_line(aes(y=run7, colour = "run 7")) +
  geom_line(aes(y=run8, colour = "run 8")) +
  
  geom_line(aes(y=stochastic_run0, colour = "stochastic run 0")) +
  geom_line(aes(y=stochastic_run1, colour = "stochastic run 1")) +
  geom_line(aes(y=stochastic_run2, colour = "stochastic run 2")) +
  geom_line(aes(y=stochastic_run3, colour = "stochastic run 3")) +
  geom_line(aes(y=stochastic_run4, colour = "stochastic run 4")) +
  geom_line(aes(y=stochastic_run5, colour = "stochastic run 5")) +
  geom_line(aes(y=stochastic_run6, colour = "stochastic run 6")) +
  geom_line(aes(y=stochastic_run7, colour = "stochastic run 7")) +
  geom_line(aes(y=stochastic_run8, colour = "stochastic run 8")) +
  geom_line(aes(y=stochastic_run9, colour = "stochastic run 9")) +
  
  geom_line(data=observed,aes(x=observed_data_matrix_x,y=observed_data_matrix_y, colour = "observed")) +
  ggtitle("volume averaged log stellar mass function") +
  xlab("log stellar mass") +
  ylab("log stellar mass function") +
  scale_y_continuous(trans='log10') +
  xlim(8.5,12)




## extract y values at x=9 -------------------------------------------------

library(pracma)
run_0_9 <- interp1(x=hist_run_0$mids,
                   y=run0_step2,
                   xi=9,
                   method="linear")
run_1_9 <- interp1(x=hist_run_0$mids,
                   y=run1_step2,
                   xi=9,
                   method="linear")
run_2_9 <- interp1(x=hist_run_0$mids,
                   y=run2_step2,
                   xi=9,
                   method="linear")
run_3_9 <- interp1(x=hist_run_0$mids,
                   y=run3_step2,
                   xi=9,
                   method="linear")
run_4_9 <- interp1(x=hist_run_0$mids,
                   y=run4_step2,
                   xi=9,
                   method="linear")
run_5_9 <- interp1(x=hist_run_0$mids,
                   y=run5_step2,
                   xi=9,
                   method="linear")
run_6_9 <- interp1(x=hist_run_0$mids,
                   y=run6_step2,
                   xi=9,
                   method="linear")
run_7_9 <- interp1(x=hist_run_0$mids,
                   y=run7_step2,
                   xi=9,
                   method="linear")
run_8_9 <- interp1(x=hist_run_0$mids,
                   y=run8_step2,
                   xi=9,
                   method="linear")


extracted_at_9 <- c(run_0_9,
                    run_1_9,
                    run_2_9,
                    run_3_9,
                    run_4_9,
                    run_5_9,
                    run_6_9,
                    run_7_9,
                    run_8_9)


## extract y values at x=10.5 ----------------------------------------------
run_0_10.5 <- interp1(x=hist_run_0$mids,
                      y=run0_step2,
                      xi=10.5,
                      method="linear")
run_1_10.5 <- interp1(x=hist_run_0$mids,
                      y=run1_step2,
                      xi=10.5,
                      method="linear")
run_2_10.5 <- interp1(x=hist_run_0$mids,
                      y=run2_step2,
                      xi=10.5,
                      method="linear")
run_3_10.5 <- interp1(x=hist_run_0$mids,
                      y=run3_step2,
                      xi=10.5,
                      method="linear")
run_4_10.5 <- interp1(x=hist_run_0$mids,
                      y=run4_step2,
                      xi=10.5,
                      method="linear")
run_5_10.5 <- interp1(x=hist_run_0$mids,
                      y=run5_step2,
                      xi=10.5,
                      method="linear")
run_6_10.5 <- interp1(x=hist_run_0$mids,
                      y=run6_step2,
                      xi=10.5,
                      method="linear")
run_7_10.5 <- interp1(x=hist_run_0$mids,
                      y=run7_step2,
                      xi=10.5,
                      method="linear")
run_8_10.5 <- interp1(x=hist_run_0$mids,
                      y=run8_step2,
                      xi=10.5,
                      method="linear")


extracted_at_10.5 <- c(run_0_10.5,
                       run_1_10.5,
                       run_2_10.5,
                       run_3_10.5,
                       run_4_10.5,
                       run_5_10.5,
                       run_6_10.5,
                       run_7_10.5,
                       run_8_10.5)




## extract stochastic y values at 9 ----------------------------------------
stochastic_run_0_9 <- interp1(x=hist_run_0$mids,
                              y=stochastic_run0_step2,
                              xi=9,
                              method="linear")
stochastic_run_1_9 <- interp1(x=hist_run_0$mids,
                              y=stochastic_run1_step2,
                              xi=9,
                              method="linear")
stochastic_run_2_9 <- interp1(x=hist_run_0$mids,
                              y=stochastic_run2_step2,
                              xi=9,
                              method="linear")
stochastic_run_3_9 <- interp1(x=hist_run_0$mids,
                              y=stochastic_run3_step2,
                              xi=9,
                              method="linear")
stochastic_run_4_9 <- interp1(x=hist_run_0$mids,
                              y=stochastic_run4_step2,
                              xi=9,
                              method="linear")
stochastic_run_5_9 <- interp1(x=hist_run_0$mids,
                              y=stochastic_run5_step2,
                              xi=9,
                              method="linear")
stochastic_run_6_9 <- interp1(x=hist_run_0$mids,
                              y=stochastic_run6_step2,
                              xi=9,
                              method="linear")
stochastic_run_7_9 <- interp1(x=hist_run_0$mids,
                              y=stochastic_run7_step2,
                              xi=9,
                              method="linear")
stochastic_run_8_9 <- interp1(x=hist_run_0$mids,
                              y=stochastic_run8_step2,
                              xi=9,
                              method="linear")
stochastic_run_9_9 <- interp1(x=hist_run_0$mids,
                              y=stochastic_run9_step2,
                              xi=9,
                              method="linear")


stochastic_extracted_at_9 <- c(stochastic_run_0_9,
                               stochastic_run_1_9,
                               stochastic_run_2_9,
                               stochastic_run_3_9,
                               stochastic_run_4_9,
                               stochastic_run_5_9,
                               stochastic_run_6_9,
                               stochastic_run_7_9,
                               stochastic_run_8_9,
                               stochastic_run_9_9)


## extract stochastic y values at 10.5 -------------------------------------
stochastic_run_0_10.5 <- interp1(x=hist_run_0$mids,
                                 y=stochastic_run0_step2,
                                 xi=10.5,
                                 method="linear")
stochastic_run_1_10.5 <- interp1(x=hist_run_0$mids,
                                 y=stochastic_run1_step2,
                                 xi=10.5,
                                 method="linear")
stochastic_run_2_10.5 <- interp1(x=hist_run_0$mids,
                                 y=stochastic_run2_step2,
                                 xi=10.5,
                                 method="linear")
stochastic_run_3_10.5 <- interp1(x=hist_run_0$mids,
                                 y=stochastic_run3_step2,
                                 xi=10.5,
                                 method="linear")
stochastic_run_4_10.5 <- interp1(x=hist_run_0$mids,
                                 y=stochastic_run4_step2,
                                 xi=10.5,
                                 method="linear")
stochastic_run_5_10.5 <- interp1(x=hist_run_0$mids,
                                 y=stochastic_run5_step2,
                                 xi=10.5,
                                 method="linear")
stochastic_run_6_10.5 <- interp1(x=hist_run_0$mids,
                                 y=stochastic_run6_step2,
                                 xi=10.5,
                                 method="linear")
stochastic_run_7_10.5 <- interp1(x=hist_run_0$mids,
                                 y=stochastic_run7_step2,
                                 xi=10.5,
                                 method="linear")
stochastic_run_8_10.5 <- interp1(x=hist_run_0$mids,
                                 y=stochastic_run8_step2,
                                 xi=10.5,
                                 method="linear")
stochastic_run_9_10.5 <- interp1(x=hist_run_0$mids,
                                 y=stochastic_run9_step2,
                                 xi=10.5,
                                 method="linear")

stochastic_extracted_at_10.5 <- c(stochastic_run_0_10.5,
                                  stochastic_run_1_10.5,
                                  stochastic_run_2_10.5,
                                  stochastic_run_3_10.5,
                                  stochastic_run_4_10.5,
                                  stochastic_run_5_10.5,
                                  stochastic_run_6_10.5,
                                  stochastic_run_7_10.5,
                                  stochastic_run_8_10.5,
                                  stochastic_run_9_10.5)





# Figure 26 ---------------------------------------------------------------

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






## emulator at 9 -----------------------------------------------------------------
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
par(mfrow=c(1,1))
plot(x,y)




xD <- matrix(data=c(x,y),ncol=2,byrow = FALSE)
D <- extracted_at_9

#which   ?*sd(D) = sd(repetitions Shaun has sent) 
#?=sd(repetitions Shaun has sent)/sd(D)
#nugget is ?^2
sd(stochastic_extracted_at_9)/sd(D) #0.05430212

#png(file="9_emulator_%d.png",width=1400,height=1200,res=200)
nugget_argument_efficient_nugget_all_graphs(xD_matrix = xD,
                                            D_vector = D,
                                            grid_length = 50,
                                            theta_val = 0.5, 
                                            sigma_val = sd(D), 
                                            E_f_val = mean(D),
                                            
                                            nugget_val = 0.05430212^2) #must input it as squared as dealing with sigma squared
#dev.off()


## emulator at 10.5 --------------------------------------------------------------
xD <- matrix(data=c(x,y),ncol=2,byrow = FALSE)
D <- extracted_at_10.5

#which   ?*sd(D) = sd(repetitions Shaun has sent) 
#?=sd(repetitions Shaun has sent)/sd(D)
#nugget is ?^2
sd(stochastic_extracted_at_10.5)/sd(D) #0.03088729

#png(file="10.5_emulator_%d.png",width=1400,height=1200,res=200)
nugget_argument_efficient_nugget_all_graphs(xD_matrix = xD,
                                            D_vector = D,
                                            grid_length = 50,
                                            theta_val = 0.5, 
                                            sigma_val = sd(D), 
                                            E_f_val = mean(D),
                                            
                                            nugget_val = 0.03088729^2) #must input it as squared as dealing with sigma squared
#dev.off()




## History Matching at 9, 10.5, max ----------------------------------------


# And again we will need our contour plotting function from before, but now we 
# adapt it so we can draw some extra contour lines using the cont_levs_lines 
# argument below, so we can highlight the boundary of the non-implausible 
# region. 
# Have a look at the following function emul_fill_cont_V2():

##############################################################################################.
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
### end define filled contour plot function for emulator output ###
##############################################################################################.

# Now input the usual colour schemes:

### Define colour schemes for standard plots ###
exp_cols <- magma
var_cols <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
diag_cols <- turbo







par(mar=c(4,4,3,3))




### history matching 9 ------------------------------------------------------
par(mar=c(4,4,3,3))


#set up

#Store the extra information we need for the History Match: 
#     the observation data 𝑧, 
#     observation errors 𝜎_𝑒= √ (  Var[𝑒]  ) and 
#     model discrepancy 𝜎_𝜖= √ (  Var[𝜖]  ), 
#where we use the standard deviation form for ease of interpretability.


### Defined Extra Objects for HM and Implausibility ###
observation9 <- interp1(x=observed_data_matrix[,"mid point"],
                        y=10^(-3)*observed_data_matrix[,"phi"],
                        xi=9,
                        method="linear")
#0.01225

z <- 0.01225 #observation recorded in real life
sigma_e <- 0.001225 #observation error
sigma_epsilon <- 0.001225 #model discrepancy

#first attempt:
#I will try out observation error 5% of value
#model discrepancy 10% of value
0.01225*0.05 #0.0006125
#initially tried 0.0006125 and 0.001225 now just trying bigger errors to get degenerate strip


#Wave 1
#Set up the 2D emulator as usual, defining the function𝑓(𝑥), the prediction 
#points 𝑥_𝑃, the design points 𝑥_𝐷 and then performing the runs 𝐷.
#Evaluate our emulator at all the 2500 input points in 𝑥_𝑃 as usual.
#Calculate Implausibility Measure Over All 2500 input points in xP
#MAYBE DELETENote that colours green and blue imply the implausibility is less than the 
# cutoff of 𝑐=3 and so these parts of the input space are deemed 
# “non-implausible”. 
# Areas where the colours are orange/red/dark red imply that 𝐼(𝑥) > 3 and 
# hence are deemed implausible and discarded from the search.



### Define wave 1 run locations: set random seed as above ###
xD_w1 <- matrix(data=c(x,y),ncol=2,byrow = FALSE)



### Define 50x50 grid of prediction points xP for emulator evaluation ###
x_grid <- seq(-0.001,1.001,len=50)
xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))



# We now fix the wave 1 runs, evaluate them and emulate as usual.
#ie wave 1 is just our starting emulator points (where we would've stopped in
#previous practicals)


### Define current run locations ###
xD <- xD_w1

### The evaluations ###
D <- extracted_at_9


### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=0.5,sigma=sd(D),E_f=mean(D),nugget=0.05430212^2)  
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

### Calculate Implausibility Measure Over All 50x50 = 2500 input points in xP ###
Imp_mat <- sqrt( (E_D_fx_mat - z)^2 / (Var_D_fx_mat + sigma_e^2 + sigma_epsilon^2) ) 

### Define colours and levels for implausibility plots ###
imp_cols <- function(n) turbo(n,begin=0.15,end=1)
imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),20,40,61)

### plot wave 1 implausibility and wave 1 runs only ###


# I use this graph in Section 10.1 Figure 33a
# See below for the graphs I use for Figure 26

#png(file="9_hm_%d.png",width=1400,height=1200,res=200)
emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3,xD=xD,x_grid=x_grid,
                  xD_col="purple",color.palette=imp_cols,main="Implausibility I(x)")
#dev.off()

#note the black lines around regions of non-implausibility!






### history matching 10.5 ---------------------------------------------------

observation10.5 <- interp1(x=observed_data_matrix[,"mid point"],
                           y=10^(-3)*observed_data_matrix[,"phi"],
                           xi=10.5,
                           method="linear")
#0.00512

z_2 <- 0.00512
sigma_e_2 <- 0.000512
sigma_epsilon_2 <- 0.000512


### Define wave 1 run locations: set random seed as above ###
xD_w1 <- matrix(data=c(x,y),ncol=2,byrow = FALSE)



### Define 50x50 grid of prediction points xP for emulator evaluation ###
x_grid <- seq(-0.001,1.001,len=50)
xP <- as.matrix(expand.grid("x1"=x_grid,"x2"=x_grid))



# We now fix the wave 1 runs, evaluate them and emulate as usual.
#ie wave 1 is just our starting emulator points (where we would've stopped in
#previous practicals)


### Define current run locations ###
xD <- xD_w1

### The evaluations ###
D <- extracted_at_10.5


### Evaluate emulator over 50x50=2500 prediction points xP and store as matrices ###
em_out <- nugget_argument_efficient_nugget_simple_BL_emulator(xP=xP,xD=xD,D=D,theta=0.5,sigma=sd(D),E_f=mean(D),nugget=0.03088729^2)  
E_D_fx_mat <- matrix(em_out[,"ExpD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 
Var_D_fx_mat <- matrix(em_out[,"VarD_f(x)"],nrow=length(x_grid),ncol=length(x_grid)) 

### Calculate Implausibility Measure Over All 50x50 = 2500 input points in xP ###
Imp_mat_2 <- sqrt( (E_D_fx_mat - z_2)^2 / (Var_D_fx_mat + sigma_e_2^2 + sigma_epsilon_2^2) )

### Define colours and levels for implausibility plots ###
imp_cols <- function(n) turbo(n,begin=0.15,end=1)
imp_levs <- c(0,seq(1,2.75,0.25),seq(3,18,2),20,40,61)




# We can now plot the implausibilities for the first output 𝐼_1(𝑥), the 
# second output 𝐼_2(𝑥) and then calculate and plot the maximised 
# implausibility defined as:

#     𝐼_𝑀(𝑥) = max_{over𝑖}    𝐼_𝑖(𝑥)

# To calculate this we collapse the matrices of implausibilities back to 
# vectors using c(), column bind them together using cbind() then maximise 
# using apply().




# Aside: 
# Note that we could have used the abind() function from the abind package here 
# which allows you to bind matrices and arrays together to produce larger 
# arrays (these are like matrices but with higher than 2 dimensions!). 
# This then makes maximising across two matrices trivial, again using apply(). 
# End Aside.





### plot implausibility for each output after wave 2 ###
#png(file="9_hm_%d.png",width=1400,height=1200,res=200)
emul_fill_cont_V2(cont_mat=Imp_mat,cont_levs=imp_levs,cont_levs_lines=3,xD=xD,
                  x_grid=x_grid,xD_col="purple",color.palette=imp_cols,
                  main="Output 9: Implausibility I_1(x)")
#dev.off()





#png(file="10.5_hm_%d.png",width=1400,height=1200,res=200)
emul_fill_cont_V2(cont_mat=Imp_mat_2,cont_levs=imp_levs,cont_levs_lines=3,xD=xD,
                  x_grid=x_grid,xD_col="purple",color.palette=imp_cols,
                  main="Output 10.5: Implausibility I_2(x)")
#dev.off()






### combine two implausibilities from f and f_2 into I_M = max_i I_i(x) ###
### convert matrices to vectors using c(), column bind using cbind(), then max using apply():
max_Implaus <- apply(cbind(c(Imp_mat),c(Imp_mat_2)),1,max)  # vector of imps. Could use abind!
#this makes matrix with 2 columns: col 1 Imp_mat vector; col 2 Imp_mat_2 vector
#then applies max to each row of matrix
#this gives vector of length 2500 (50x50 grid so 2500 pts)
max_Imp_mat <- matrix(max_Implaus,nrow=length(x_grid),ncol=length(x_grid)) # convert to matrix
#converts 2500 length vector of max implausibilities to 50x50 matrix


#png(file="max_hm_%d.png",width=1400,height=1200,res=200)
emul_fill_cont_V2(cont_mat=max_Imp_mat,cont_levs=imp_levs,cont_levs_lines=3,
                  xD=xD,x_grid=x_grid,xD_col="purple",color.palette=imp_cols,
                  main="Max Implausibility I_M(x)")
#dev.off()










# Figure 27 ---------------------------------------------------------------


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



## The following produce the best design for given batch size --------------
# They are commented out so it is easier to uncomment batch size interested in and ctrl A run

#I will use a ceiling of 5000 this time...


### 4 -----------------------------------------------------------------------
# png(file="pts4_5000_%d.png",width=1400,height=1200,res=200)
# hm_pts4_5000 <- hm_everything_function_sigmoidal(want_all_var_graphs=FALSE,
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
#                                                   batch_size=4, #4, 6, 8 or 10
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
#                                                   ceiling_value=5000, #NOT LOGGED max total run time cannot exceed
# 
# 
#                                                   #NEW:
#                                                   z_observation=z_10.5,
#                                                   sigma_e=sigma_e_10.5,
#                                                   sigma_epsilon=sigma_epsilon_10.5
# )
# dev.off()
# saveRDS(hm_pts4_5000, file="hm_pts4_5000.RData")

### 5 -----------------------------------------------------------------------
# png(file="pts5_5000_%d.png",width=1400,height=1200,res=200)
# hm_pts5_5000 <- hm_everything_function_sigmoidal(want_all_var_graphs=FALSE,
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
#                                                   batch_size=5, #4, 6, 8 or 10
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
#                                                   ceiling_value=5000, #NOT LOGGED max total run time cannot exceed
# 
# 
#                                                   #NEW:
#                                                   z_observation=z_10.5,
#                                                   sigma_e=sigma_e_10.5,
#                                                   sigma_epsilon=sigma_epsilon_10.5
# )
# dev.off()
# saveRDS(hm_pts5_5000, file="hm_pts5_5000.RData")


### 6 -----------------------------------------------------------------------
# png(file="pts6_5000_%d.png",width=1400,height=1200,res=200)
# hm_pts6_5000 <- hm_everything_function_sigmoidal(want_all_var_graphs=FALSE,
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
#                                                   ceiling_value=5000, #NOT LOGGED max total run time cannot exceed
# 
# 
#                                                   #NEW:
#                                                   z_observation=z_10.5,
#                                                   sigma_e=sigma_e_10.5,
#                                                   sigma_epsilon=sigma_epsilon_10.5
# )
# dev.off()
# saveRDS(hm_pts6_5000, file="hm_pts6_5000.RData")



### 7 -----------------------------------------------------------------------
# png(file="pts7_5000_%d.png",width=1400,height=1200,res=200)
# hm_pts7_5000 <- hm_everything_function_sigmoidal(want_all_var_graphs=FALSE,
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
#                                                   batch_size=7, #4, 6, 8 or 10
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
#                                                   ceiling_value=5000, #NOT LOGGED max total run time cannot exceed
# 
# 
#                                                   #NEW:
#                                                   z_observation=z_10.5,
#                                                   sigma_e=sigma_e_10.5,
#                                                   sigma_epsilon=sigma_epsilon_10.5
# )
# dev.off()
# saveRDS(hm_pts7_5000, file="hm_pts7_5000.RData")



### 8 -----------------------------------------------------------------------
# png(file="pts8_5000_%d.png",width=1400,height=1200,res=200)
# hm_pts8_5000 <- hm_everything_function_sigmoidal(want_all_var_graphs=FALSE,
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
#                                                   batch_size=8, #4, 6, 8 or 10
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
#                                                   ceiling_value=5000, #NOT LOGGED max total run time cannot exceed
# 
# 
#                                                   #NEW:
#                                                   z_observation=z_10.5,
#                                                   sigma_e=sigma_e_10.5,
#                                                   sigma_epsilon=sigma_epsilon_10.5
# )
# dev.off()
# saveRDS(hm_pts8_5000, file="hm_pts8_5000.RData")



### 9 -----------------------------------------------------------------------
# png(file="pts9_5000_%d.png",width=1400,height=1200,res=200)
# hm_pts9_5000 <- hm_everything_function_sigmoidal(want_all_var_graphs=FALSE,
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
#                                                   batch_size=9, #4, 6, 8 or 10
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
#                                                   ceiling_value=5000, #NOT LOGGED max total run time cannot exceed
# 
# 
#                                                   #NEW:
#                                                   z_observation=z_10.5,
#                                                   sigma_e=sigma_e_10.5,
#                                                   sigma_epsilon=sigma_epsilon_10.5
# )
# dev.off()
# saveRDS(hm_pts9_5000, file="hm_pts9_5000.RData")









# Figure 28 ---------------------------------------------------------------

## now analysing the results... --------------------------------------------

#check Tc and v90w2 for them all
#4
pts4 <- readRDS("hm_pts4_5000.RData")
pts4$`Best Design v90 W2` #1.243572e-06
pts4$`Best Design Tc` #3315.208
#5
pts5 <- readRDS("hm_pts5_5000.RData")
pts5$`Best Design v90 W2` #7.587337e-07
pts5$`Best Design Tc` #3882.15
#6
pts6 <- readRDS("hm_pts6_5000.RData")
pts6$`Best Design v90 W2` #5.73303e-07
pts6$`Best Design Tc` #4534.503
#7
pts7 <- readRDS("hm_pts7_5000.RData")
pts7$`Best Design v90 W2` #8.380537e-07
pts7$`Best Design Tc` #4991.317
#8
pts8 <- readRDS("hm_pts8_5000.RData")
pts8$`Best Design v90 W2` #4.183402e-06
pts8$`Best Design Tc` #4848.324
#9
pts9 <- readRDS("hm_pts9_5000.RData")
pts9$`Best Design v90 W2` #4.475702e-06
pts9$`Best Design Tc` #4996.593



## best design -------------------------------------------------------------

# Did not use this one in report - see the graph below
#png(file="hm_optimal_number_of_points.png",width=1200,height=1315,res=200)

plot(x=c(4,5,6,7,8,9),y=c(1.243572e-06,
                          7.587337e-07,
                          5.73303e-07,
                          8.380537e-07,
                          4.183402e-06,
                          4.475702e-06),
     xlab="Number of Points",
     ylab="V90",
     main="Finding the Optimal Number of Points",
     pch = 19
)
lines(x=c(4,5,6,7,8,9),y=c(1.243572e-06,
                           7.587337e-07,
                           5.73303e-07,
                           8.380537e-07,
                           4.183402e-06,
                           4.475702e-06))

#dev.off()


#new graph showing 0 on y-axis - used this in report
#png(file="new0_hm_optimal_number_of_points.png",width=1200,height=1315,res=200)
plot(x=c(4,5,6,7,8,9),y=c(1.243572e-06,
                          7.587337e-07,
                          5.73303e-07,
                          8.380537e-07,
                          4.183402e-06,
                          4.475702e-06),
     xlab="Number of Points",
     ylab="V90",
     main="Finding the Optimal Number of Points",
     pch = 19,
     ylim = c(0,4.475702e-06)
)
lines(x=c(4,5,6,7,8,9),y=c(1.243572e-06,
                           7.587337e-07,
                           5.73303e-07,
                           8.380537e-07,
                           4.183402e-06,
                           4.475702e-06))
#dev.off()





# End of Section 8 --------------------------------------------------------
