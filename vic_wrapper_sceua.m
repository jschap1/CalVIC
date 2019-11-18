% VIC wrapper for the shuffled complex algorithm
% October 31, 2019 JRS
%
% INPUTS
% x = parameters
% icall = number of the current model run
%
% OUTPUTS
% f = error metric (NSE)

function f = vic_wrapper_sceua(nopt, x, control_params, icall)

icall = icall + 1; % using 1-indexing instead of 0-indexing

% Create directory for VIC model outputs
control_params.vic_out_dir = [control_params.vic_out_dir num2str(icall)];
mkdir(control_params.vic_out_dir)
disp(['Created directory ', control_params.vic_out_dir])

% Create directory for glacier melt model outputs
control_params.glacier_out_dir = [control_params.glacier_out_dir num2str(icall)];
mkdir(control_params.glacier_out_dir)
disp(['Created directory ', control_params.glacier_out_dir])

% Create directory for routing model inputs
control_params.rout_in_dir = [control_params.rout_in_dir num2str(icall)];
mkdir(control_params.rout_in_dir)
disp(['Created directory ', control_params.rout_in_dir])

% Create directory for routing model outputs
control_params.rout_out_dir = [control_params.rout_out_dir num2str(icall)];
mkdir(control_params.rout_out_dir)
disp(['Created directory ', control_params.rout_out_dir])

%% Other inputs (not used)

% Template global parameter file
% filenames.global_param_template = '/Volumes/HD3/SWOTDA/Calibration/Pandoh/par_2015_WB.txt';
% filenames.global_param_file = '/Volumes/HD3/SWOTDA/Calibration/onecell/par_2015_WB.txt';
% Initial routing model parameter file
% filenames.init_rout_pars = '/Volumes/HD3/SWOTDA/Calibration/Pandoh/rout_param.txt';
% Initial soil parameter file (first guess)
% filenames.soil_initial = '/Volumes/HD3/SWOTDA/Calibration/Pandoh/soils_pandoh.txt';
% filenames.init_soil_pars = '/Volumes/HD3/SWOTDA/Calibration/onecell/soils_onecell_init.txt';
% Discharge observations
% filenames.discharge_obs = '/Volumes/HD3/SWOTDA/Calibration/Pandoh/pandoh_inflow.txt';
% Glacier fraction raster that matches the basin mask
% filenames.glacier_fract_maps
% The outputs and soil parameter files are saved separately for each model run
% Location to save the routing model inputs
% Location to save the routing model outputs
% filenames.soil_param = filenames.soil_initial;

control_params.soil_param = control_params.init_soil_pars;

%% Write soil parameter file using the parameters x

% Load soil parameters
if icall==1
     soils = load(control_params.init_soil_pars);
%      control_params.soil_param = control_params.init_soil_pars;
else
    soils = load(control_params.soil_param);
  
  % Replace the old parameters with the new parameters
%   soils(:,25) = x;
%   soils(:,25) = mean(x);

    soils(:,5) = x(1); % b_infilt
    soils(:,6) = x(2); % ds
    soils(:,8) = x(3); % ws

%     soils(:,5) = x(1); % b_infilt
%     soils(:,6) = x(2); % ds
%     soils(:,7) = x(3); % dsmax
%     soils(:,8) = x(4); % ws
    
%     soils(:,24) = x(5); % layer 2 soil depth
%     soils(:,25) = x(6); % layer 3 soil depth
        
%     soils(:,5) = x;
    
    % Write out the updated soil parameter file
    soil_path = fileparts(control_params.soil_param);
    control_params.soil_param = [soil_path '/soils_' num2str(icall) '.txt'];
    write_soils(5, soils, control_params.soil_param, '3l');
  
end


%% Write global parameter file

% Read the template global parameter file
% filenames.global_param = filenames.global_param_file;
A = read_global_param_file(control_params.global_param_file);

% if icall==1
%   A = read_global_param_file(filenames.global_param_template);
% else
%   filenames.global_param = ['/Volumes/HD3/SWOTDA/Calibration/Pandoh/global_param_' num2str(icall-1) '.txt'];
%   A = read_global_param_file(filenames.global_param);
% end

% Change the reference to the soil parameter file
A{153} = ['SOIL ' control_params.soil_param];

% Change output directory
A{181} = ['LOG_DIR ' control_params.vic_out_dir '/'];
A{182} = ['RESULT_DIR ' control_params.vic_out_dir '/'];

% Write global parameter file
% filenames.global_param = ['/Volumes/HD3/SWOTDA/Calibration/Pandoh/global_param_' num2str(icall) '.txt'];
write_global_param_file(A, control_params.global_param_file)

%% Run VIC

% Check if outputs already exist to avoid running VIC more than necessary
fn1 = dir(fullfile(control_params.vic_out_dir, '*.txt'));
if length(fn1)>0
    
    disp('VIC has already been run for this iteration')
    
else
    
    % Split the VIC run over multiple processors (for Hoffman2 UGE only)
    if control_params.n_proc > 1

        vic_exec_file = set_up_parallel(control_params);
        vic_run_command = ['qsub -cwd -V -N PAN -l h_data=1024M,h_rt=00:30:00 -M $HOME -m bea -t 1-' num2str(control_params.n_proc) ':1 ' vic_exec_file];
        system(vic_run_command)

%         qsub -cwd -V -N PJ -l h_data=1024M,h_rt=01:00:00 -M $HOME -m bea -t 1-100:1 myFuncWrapper.sh
        
    else    
    
        vic_run_command = [control_params.vic_command ' -g ' control_params.global_param_file];
        system(vic_run_command)
        
    end
    
end

%% Adjust runoff for glacier melt contribution and write routing model input files

ddf = 7; % degree-day factor (mm/K/day)
plotflag = 0;
add_glacier_contribution(ddf, control_params.glacier_fract_map, ...
    control_params.vic_out_dir, control_params.rout_in_dir, ...
    control_params.glacier_out_dir, plotflag)

%% Routing model

% Write routing model input file
B = read_global_param_file(control_params.rout_param_file);
B{19} = fullfile(control_params.rout_in_dir, 'fluxes_');
B{22} = [control_params.rout_out_dir '/'];
write_global_param_file(B, control_params.rout_param_file)

% Run routing model
rout_run_command = [control_params.rout_command ' ' control_params.rout_param_file];
system(rout_run_command)

% Remove existing unit hydrograph files
!rm *.uh_s

%% Compute error metric

% Load gage discharge
obs = dlmread(control_params.discharge_obs, '\t', 1, 0);
% obs = dlmread(control_params.discharge_obs, '\t', 3, 0);
Qobs = obs(:,4);
time_obs = datetime(obs(:,1), obs(:,2), obs(:,3));
Qobs(Qobs<0) = NaN; % remove negative "observed" values of discharge

% obs = dlmread(filenames.obs_discharge, '\t', 1, 0);
% Qobs = obs(:,4)/1000; % cfs (dividing by 1000 for the one-cell test)
% Qobs = obs(:,4) + obs(:,5); % combining runoff and baseflow for the one-cell test
% Qobs = Qobs + 10; % add a bias, just because

% Load discharge calculated with routing model
% rout = dlmread(['/Volumes/HD3/SWOTDA/Calibration/Pandoh/vic/out_' num2str(icall) '/fluxes_31.84375_77.03125.txt'], '\t', 3, 0);
% rout = dlmread(['/Volumes/HD3/SWOTDA/Calibration/onecell/vic_out_' num2str(icall) '/fluxes_31.84375_77.03125.txt'], '\t', 3, 0);
% Q = rout(:,4) + rout(:,5);

rout = load(fullfile(control_params.rout_out_dir, 'PAN1 .day'));
time_rout = datetime(rout(:,1), rout(:,2), rout(:,3));
Q = rout(:,4);

% Make sure the two time series are consistent
start_time = max(time_rout(1), time_obs(1));
end_time = min(time_rout(end), time_obs(end));
Qobs = Qobs(find(time_obs == start_time):find(time_obs == end_time));
Q = Q(find(time_rout == start_time):find(time_rout == end_time));
time_merged = time_rout(find(time_rout == start_time):find(time_rout == end_time));

plotflag = 0;
if plotflag
    figure
    jsplot(time_merged, Q, 'Flow', 'Time', 'Discharge (cfs)', 18)
    hold on
    plot(time_merged, Qobs)
    legend('Predicted','Observed','Location','NW')
    saveas(gcf, fullfile(control_params.rout_out_dir, ['discharge_plot_iter_' num2str(icall), '.png']))
end

if length(Q) ~= length(Qobs)
  error('Lengths of discharge records are not equal')
end

%% Calculate error metric, ignore the first few time steps (spin-up)

Qobs = Qobs(control_params.n_spinup:end);
Q = Q(control_params.n_spinup:end);

error_metric = 'RMSE';
if strcmp(error_metric, 'NSE')
  f = myNSE(Qobs, Q);
  disp(['NSE for iteration ' num2str(icall) ' is ' num2str(f)])
elseif strcmp(error_metric, 'RMSE')
  f = myRMSE(Qobs, Q);
  disp(['RMSE for iteration ' num2str(icall) ' is ' num2str(f)])
end

return
