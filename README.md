# bayesCONUS
 The computer code used to run the Bayesian model and produce the results in Piecuch (2025), The Rate of U.S. Coastal Sea-Level Rise Doubled in the Past Century, written in the MATLAB software environment 

bayesCONUS

README file last updated by Christopher Piecuch, cpiecuch-at-whoi-dot-edu, Sun Oct 19 2025

Basic description

Citation

This code was generated to produce the results presented in the main text of:

Piecuch, C. G.,  “The Rate of U.S. Coastal Sea-Level Rise Doubled in the Past Century”, AGU Advances, 2025AV002018.

Please cite this reference when using this code. 

Contents

Text files
•	Copyright: copyright statement
•	License: license statement

MATLAB .m files
•	bayes_run_analyses.m: this is the main driver code of the model.  Simply execute “bayes_run_analyses” in the MATLAB Command Window from this directory, and this code should run “out of the box” and produce the figures in the paper after a few hours.
•	bayes_main_code.m: This performs the Bayesian inference at and based on the tide-gauge data
•	bayes_posterior_prediction.m: This performs posterior prediction, taking the solutions from bayes_main_code.m at tide-gauge locations and makes inference at the regular 0.5-degree CONUS grid
•	bayes_make_figures.m: This makes the figures shown in the paper.
•	bayes_plot_residuals.m: This plots residual differences between Bayesian model solution and tide-gauge data.
•	delete_burn_in.m: delete “burn-in” (or “warm-up”) transients from model solution
•	determine_clusters.m: create matrix used in spatial covariance structure of relative sea level fluctuations
•	EarthDistances.m: compute distances between latitude and longitude points on a spherical Earth
•	grab_coastal_cells.m: create regular grid of coastal locations
•	init_vals_pickup.m: initialize output from pickup file
•	initialize_output.m: initialize output
•	prepare_tgr_data.m: format tide gauge data and bring into Bayesian model workspace
•	randraw.m: draw random value from any of a number of distributions
•	set_hyperparameters.m: set hyperparameter values (i.e., parameters of the prior distributions)
•	set_initial_values.m: set initial values of model solutions
•	update_all_arrays.m: update model solutions

Subdirectories
•	rlr_annual_20250802: tide gauge records of relative sea level

Note on usage of code and MATLAB release
This code was written and executed using the MATLAB release version MATLAB_R2024a. 
