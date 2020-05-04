CalVIC: A Matlab toolbox for calibrating the Variable Infiltration Capacity (VIC) model

run_sceua.m
	reads parameter file
	calls sceua.m
sceua.m
	sets up complexes
	calls vic_wrapper_sceua.m at each sample point
	given the rmse value returned by vic_wrapper_sceua, calls cceua.m

cceua.m 
	evolves each complex to produce new guesses
	shuffles points among complexes	

vic_wrapper_sceua.m
	runs the vic model, distributing the job across Hoffman2 nodes
	adds glacier contribution to runoff
	runs the routing model to predict streamflow
	uses predicted streamflow to calculate function value (f = RMSE)
	
	