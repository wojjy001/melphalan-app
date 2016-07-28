##Some explanation for the plots that are in this folder##
Plots titled:
	- R_PDSIM_OrigData_100.png
	- R_PKSIM_OrigData_100.png
			- These are plots created by an R script that mimics a Visual Predictive Check by
			simulating your index dataset 100 times
			- Blue circles are the observed concentrations and counts
			- Observed data are also summarised by the median (red solid line) and 5th and 95th
			percentiles (red dashed lines)
			- Simulated data are summarised by the 95% prediction intervals for the median (red 
			shaded ribbon) and 95% prediction intervals for the 5th and 95th percentiles (blue shaded 
			ribbons)
			- As individuals were sampled at different time points, observations were binned into 10 		
			groups for calculation of summary statistics
	
	- NONMEM_PDSIM_OrigData_100.png
	- NONMEM_PKSIM_OrigData_100.png
			- These plots represent 100 simulations of the index dataset by NONMEM (Visual Predictive 			
			Check style as the previous R simulations)
			* Plots produced by the R script and NONMEM are in good agreement


Plots titled:
	- R_PDSIM_1000.png
	- R_PKSIM_1000.png
			- These are plots created by the "Melphalan_Simulation_Model.R" script
			- 1000 individuals with the same characteristics were simulated
			- Melphalan concentrations and neutrophil counts were calculated based on the simulated
			individual parameters and differential equations
			- Plots show the median (solid line) and 5th and 95th percentiles (shaded ribbon) for
			predictions
			
	- NONMEM_PDSIM_1000.png
	- NONMEM_PKSIM_1000.png
			- When a patient with the same characteristics as the previous R simulations is also
			simulated 1000 times using NONMEM, the resultant plots are given
			* Plots produced by the "Melphalan_Simulation_Model.R" script and NONMEM are in good
			agreement			
			
			
			
			
			
			
			
			
			
			
			
			