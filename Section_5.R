
# Section 5: Cosmological Simulations -------------------------------------

# Figures 18 and 19

# Functions ---------------------------------------------------------------

# Figure 18 ---------------------------------------------------------------
#run 4 is central 
run_4_dm_map <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/Run_4/density_maps/dm_map.txt")
matrix_run_4_dm_map <- as.matrix(run_4_dm_map)
par(mfrow=c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1))

#png(file="central_point_dark_matter.png",width=1200,height=1315,res=200)
image(log10(matrix_run_4_dm_map),
      yaxt="n",
      xaxt="n",
      main="Dark Matter Density Map \n Run Time: 1333 CPU Hours")
#dev.off()





# Figure 19 ---------------------------------------------------------------

run_2_gas_map <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/Run_2/density_maps/gas_map.txt")
matrix_run_2_gas_map <- as.matrix(run_2_gas_map)
par(mfrow=c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1))

#png(file="cheap_gas_map.png",width=1200,height=1315,res=200)
image(log10(matrix_run_2_gas_map),
      yaxt="n",
      xaxt="n",
      main="Gas Density Map \n Cheaper: 428 CPU Hours")
#dev.off()

run_0_gas_map <- read.table("/Users/amberthoms/Library/CloudStorage/OneDrive-DurhamUniversity/Year 4/Project/Information given/Shaun/OneDrive_1_20-01-2025/Stellar_AGN_feedback_2D_v2/Run_0/density_maps/gas_map.txt")
matrix_run_0_gas_map <- as.matrix(run_0_gas_map)
par(mfrow=c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1))

#png(file="expensive_gas_map.png",width=1200,height=1315,res=200)
image(log10(matrix_run_0_gas_map),
      yaxt="n",
      xaxt="n",
      main="Gas Density Map \n Expensive: 1278 CPU Hours")
#dev.off()



# End of Section 5 --------------------------------------------------------
