# Parameter file for CalVIC
# Line numbers of this file must not be changed (given my current, primitive Matlab parameter file reader)

# Locations of software packages ----------------------------

vic_command	/Volumes/HD3/SWOTDA/Software/VIC-VIC.5.1.0.rc1/vic/drivers/classic/vic_classic.exe
#vic_command	/u/home/j/jschaper/vic/VIC5/vic/drivers/classic/vic_classic.exe

rout_command	/Volumes/HD3/SWOTDA/Software/route_1.0/src/rout
#rout_command	/u/home/j/jschaper/route_1.0/src/rout

sce_dir	/Volumes/HD3/SWOTDA/Calibration/CalVIC/sce_matlab
#sce_dir	/u/home/j/jschaper/Codes/CalVIC/sce_matlab

# Input files -----------------------------------------------

global_param_file	/Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/par_2013-2019_EB.txt
#global_param_file	/u/home/j/jschaper/vic/Pandoh/par_2015_WB.txt

init_soil_pars	/Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/soils_pandoh_lumped.txt
#init_soil_pars	/u/home/j/jschaper/vic/Pandoh/soils_pandoh.txt

rout_param_file /Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/rout_param.txt
#rout_param_file /u/home/j/jschaper/vic/Pandoh/rout_param.txt

discharge_obs	/Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/pandoh_inflow.txt
#discharge_obs	/u/home/j/jschaper/vic/Pandoh/pandoh_inflow.txt

glacier_fract_map	/Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/pandoh_glacier_fract.mat
#glacier_fract_map	/u/home/j/jschaper/vic/Pandoh/pandoh_glacier_fract.tif


# Control parameters ----------------------------------------

vic_out_dir	/Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/out/vic_out_
#vic_out_dir	/u/home/j/jschaper/vic/Pandoh/out/vic_out_

glacier_out_dir	/Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/out/glacier_out_
#glacier_out_dir	/u/home/j/jschaper/vic/Pandoh/out/glacier_out_

rout_in_dir	/Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/out/rout_in_
#rout_in_dir	/u/home/j/jschaper/vic/Pandoh/out/rout_in_

rout_out_dir	/Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/out/rout_out_
#rout_out_dir	/u/home/j/jschaper/vic/Pandoh/out/rout_out_


# SCE-UA parameters --------------------------------------------
n_spinup	60
n_complexes	1
max_iter  1
kstop 10
pcento  1
peps  0.01
iseed -1
iniflg  1

soil_pars 5,6,7,8
par_bounds  /Volumes/HD3/SWOTDA/Calibration/CalVIC/parameter_bounds.txt

# Parallelization parameters --------------------------------------------
use_hoffman	0
jobname	MYJOB
data_per_job	1024
time_per_job	0.5
n_proc	20
wait_time	7
meta_out_dir	/Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/out/meta/

# Glacier model parameters -------------------------------------------
add_glacier_contribution	1
ddf	1

# Other parameters -------------------------------------------
objective	RMSE
all_outputs	1
lumped  1
basin_area 5280
calibrate_ddf 1
