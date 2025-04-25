# R Code to Reproduce Plots in Report
This repository contains the R code used to generate the Figures presented in my report, titled **"Incorporating Run Time into the Emulation of Expensive Computer Models, with Applications to Cosmological Simulations."**

**Author**: Amber Thoms

**Supervisor**: Professor Ian Vernon

**Collaborator**: Dr Shaun Brown

For any questions, feel free to contact me via **amber.s.thoms@durham.ac.uk**.

## Repository Structure
Each R script corresponds to a specific Section in the report and is divided into segments to reproduce each Figure.
In addition, there are two folders titled "Stellar_AGN_feedback_2D_v2" and "Stochasticity_tests".
These folders contain the run locations, run times and halo catalogues for the EAGLE cosmological simulations performed by Dr Shaun Brown.
"Stellar_AGN_feedback_2D_v2" contains information about the evaluations of the Orthogonal Latin Hypercube and Wave 2 design (and additional run for completion). 
"Stochasticity_tests" contains information about the stochastic repetitions of a central point x1=x2=0.48698 which has been rerun 10 times.

These scripts were developed using **R version 4.4.3 (2025-02-28) – "Trophy Case"**.

## Reproducing a Figure
To recreate a Figure:

1. Open the relevant .R script in RStudio.
2. Ensure any required data files are in the correct directory.
3. Run the script.
4. The output will be displayed in the plots plane.
5. To save the output (in the dimensions shown in the report) uncomment the relevant png(file="???.png",width=???,height=???,res=???) and dev.off() lines and rerun the script.
