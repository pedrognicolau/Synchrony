# Scripts for Seasonality, density dependence and spatial population synchrony
 
This project corresponds to the scripts and data files used for the manuscript "Seasonality, density dependence and spatial population synchrony" by Pedro G. Nicolau, Rolf A. Ims, Sigrunn H. Sørbye & Nigel G. Yoccoz

The folder structure is as follows.

\Sychrony
	\Data: important data files used to reproduce the code, and raw data 
		\count data: 
			\raw_chist_00_20.csv and \CR_CaptHist_KMdata_18_20.csv
				raw files of capture histories between 2000 and 2020, must be combined
		\Abundance_CRINLA_zeroc3.rds
			formatted data to be able to be read by INLA as multinomial, removing the animals only captured in day 3
			in the early years when they were still trapped. 
		     
		\gps.csv
			gps coordinates of stations
		\distance_matrix.csv
			distance matrix from gps.csv
		\winter_precipitation.rds & \weather_vars.rds & \popmodels_residuals.rds & \data4synchrony.rds
			all are files needed for the last step, all obtained in the scripts 
	
	
	\Scripts: scripts necessary for analysis, numbered by order of sequence (some with subnumbering). 
		\0_Important_functions.R
			important functions to compute Bayesian R^2 and smoothed correlograms
		\01_process_CR_abundance_stacked.R
			01 contains processing for 1
		\1_populationmodels_and_residuals.R
			Population models and residuals (Step 3)
		\2_Synchrony_estimation_and_plots.R
			Estimate and plot synchrony (Step 4)
		\03_weather_processing.R
			Process weather data necessary for 3
		\3_weather_effects.R
			Model weather synchrony vs pop synchrony (Step 5)

	Plots: diverse plots used (or not) in the manuscript; not needed for analysis
