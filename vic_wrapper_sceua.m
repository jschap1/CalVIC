% VIC wrapper for the shuffled complex algorithm
% October 31, 2019 JRS
%
% INPUTS
% x = parameters
% icall = number of the current model run
%
% OUTPUTS
% f = error metric (NSE)

function f = vic_wrapper_sceua(nopt, x, icall)

icall = icall + 1; % using 1-indexing instead of 0-indexing

%% Other inputs

% Template global parameter file
% filenames.global_param_template = '/Volumes/HD3/SWOTDA/Calibration/Pandoh/par_2015_WB.txt';
filenames.global_param_template = '/Volumes/HD3/SWOTDA/Calibration/onecell/par_2015_WB.txt';

% Initial routing model parameter file
filenames.rout_param = '/Volumes/HD3/SWOTDA/Calibration/Pandoh/rout_param.txt';

% Initial soil parameter file (first guess)
% filenames.soil_initial = '/Volumes/HD3/SWOTDA/Calibration/Pandoh/soils_pandoh.txt';
filenames.soil_initial = '/Volumes/HD3/SWOTDA/Calibration/onecell/soils_onecell_init.txt';

% Discharge observations
filenames.obs_discharge = '/Volumes/HD3/SWOTDA/Calibration/Pandoh/pandoh_inflow.txt';

% Glacier fraction raster that matches the basin mask
filenames.glacier_fract = '/Volumes/HD3/SWOTDA/FDT/v10282019/pandoh_glacier_fract.tif';

% The outputs and soil parameter files are saved separately for each model run

% Location to save the VIC model outputs
filenames.vic_results = ['/Volumes/HD3/SWOTDA/Calibration/onecell/vic_out_' num2str(icall) '/'];
% filenames.vic_results = ['/Volumes/HD3/SWOTDA/Calibration/Pandoh/vic/out_' num2str(icall) '/'];
mkdir(filenames.vic_results)
disp(['Created directory ', filenames.vic_results])

% Location to save the glacier melt model outputs
filenames.glacier_contribution = ['/Volumes/HD3/SWOTDA/Calibration/Pandoh/glacier_contribution_' num2str(icall) '.mat'];

% Location to save the routing model inputs
filenames.rout_in = ['/Volumes/HD3/SWOTDA/Calibration/Pandoh/rout/in_' num2str(icall)  '/'];
mkdir(filenames.rout_in)
disp(['Created directory ', filenames.rout_in])

% Location to save the routing model outputs
filenames.rout_results = ['/Volumes/HD3/SWOTDA/Calibration/Pandoh/rout/out_' num2str(icall) '/'];
mkdir(filenames.rout_results)
disp(['Created directory ', filenames.rout_results])

filenames.soil_param = filenames.soil_initial;

%% Write soil parameter file using the parameters x

if icall==1
  soils = load(filenames.soil_initial);
  filenames.soil_param = filenames.soil_initial;
else
  
  % Load the old soil parameters  
  soils = load(filenames.soil_param);
  
  % Replace the old parameters with the new parameters
%   soils(:,25) = x;
%   soils(:,25) = mean(x);

    soils(:,5) = x(1); % b_infilt
    soils(:,6) = x(2); % ds
    soils(:,7) = x(3); % dsmax
    soils(:,8) = x(4); % ws
    
%     soils(:,24) = x(5); % layer 2 soil depth
%     soils(:,25) = x(6); % layer 3 soil depth
        
%     soils(:,5) = x;
    
  % Write out the updated soil parameter file
%   filenames.soil_param = ['/Volumes/HD3/SWOTDA/Calibration/Pandoh/soils_pandoh_' num2str(icall) '.txt'];
  filenames.soil_param = ['/Volumes/HD3/SWOTDA/Calibration/onecell/soils_onecell_' num2str(icall) '.txt'];
  
  write_soils(5, soils, filenames.soil_param, '3l');
  
end


%% Write global parameter file

% Read the template global parameter file
filenames.global_param = filenames.global_param_template;
A = read_global_param_file(filenames.global_param);

% if icall==1
%   A = read_global_param_file(filenames.global_param_template);
% else
%   filenames.global_param = ['/Volumes/HD3/SWOTDA/Calibration/Pandoh/global_param_' num2str(icall-1) '.txt'];
%   A = read_global_param_file(filenames.global_param);
% end

% Change the reference to the soil parameter file
A{153} = ['SOIL ' filenames.soil_param];

% Change output directory
A{181} = ['LOG_DIR ' filenames.vic_results];
A{182} = ['RESULT_DIR ' filenames.vic_results];

% Write global parameter file
% filenames.global_param = ['/Volumes/HD3/SWOTDA/Calibration/Pandoh/global_param_' num2str(icall) '.txt'];
write_global_param_file(A, filenames.global_param)

%% Run VIC

% Check if outputs already exist to avoid running VIC more than necessary

fn1 = dir(fullfile(filenames.vic_results, '*.txt'));
if length(fn1)>0
    disp('VIC has already been run for this iteration')
else
    vic_run_command = ['/Volumes/HD3/SWOTDA/Software/VIC5/vic/drivers/classic/vic_classic.exe -g ' filenames.global_param];
    system(vic_run_command)
end

%% Adjust runoff for glacier melt contribution and write routing model input files

% ddf = 7; % degree-day factor (mm/K/day)
% plotflag = 0;
% 
% add_glacier_contribution(ddf, filenames.glacier_fract, filenames.vic_results, filenames.rout_in, filenames.glacier_contribution, plotflag)

%% Routing model

% % Write routing model input file
% B = read_global_param_file(filenames.rout_param);
% B{19} = [filenames.rout_in 'fluxes_'];
% B{22} = filenames.rout_results;
% write_global_param_file(B, filenames.rout_param)
% 
% % Run routing model
% rout_run_command = ['/Volumes/HD3/SWOTDA/Software/route_1.0/src/rout ' filenames.rout_param];
% system(rout_run_command)
% !rm *.uh_s

%% Compute error metric

% Load gage discharge
% obs = dlmread(filenames.obs_discharge, '\t', 1, 0);
obs = dlmread('/Volumes/HD3/SWOTDA/Calibration/onecell/true/fluxes_31.84375_77.03125.txt', '\t', 3, 0);
% obs = dlmread(filenames.obs_discharge, '\t', 1, 0);
% Qobs = obs(:,4)/1000; % cfs (dividing by 1000 for the one-cell test)
Qobs = obs(:,4) + obs(:,5); % combining runoff and baseflow for the one-cell test
Qobs(Qobs<0) = NaN; % remove negative "observed" values of discharge
% Qobs = Qobs + 10; % add a bias, just because

time_obs = datetime(obs(:,1), obs(:,2), obs(:,3));

% Load discharge calculated with routing model

% rout = dlmread(['/Volumes/HD3/SWOTDA/Calibration/Pandoh/vic/out_' num2str(icall) '/fluxes_31.84375_77.03125.txt'], '\t', 3, 0);
rout = dlmread(['/Volumes/HD3/SWOTDA/Calibration/onecell/vic_out_' num2str(icall) '/fluxes_31.84375_77.03125.txt'], '\t', 3, 0);
% rout = load(fullfile(filenames.rout_results, 'PAN1 .day'));
Q = rout(:,4) + rout(:,5);
time_rout = datetime(rout(:,1), rout(:,2), rout(:,3));

% Make sure the two time series are consistent
start_time = max(time_rout(1), time_obs(1));
end_time = min(time_rout(end), time_obs(end));
Qobs = Qobs(find(time_obs == start_time):find(time_obs == end_time));
Q = Q(find(time_rout == start_time):find(time_rout == end_time));
time_merged = time_rout(find(time_rout == start_time):find(time_rout == end_time));

plotflag = 1;
if plotflag
    figure
    jsplot(time_merged, Q, 'Flow', 'Time', 'Discharge (cfs)', 18)
    hold on
    plot(time_merged, Qobs)
    legend('Predicted','Observed','Location','NW')
    saveas(gcf, fullfile(filenames.rout_results, ['discharge_plot_iter_' num2str(icall), '.png']))
end

if length(Q) ~= length(Qobs)
  warning('Lengths are not equal, pausing')
  pause;
end

%% Calculate error metric, ignore the first few time steps (spin-up)

n_spinup = 60;
Qobs = Qobs(n_spinup:end);
Q = Q(n_spinup:end);

error_metric = 'RMSE';
if strcmp(error_metric, 'NSE')
  f = myNSE(Qobs, Q);
  disp(['NSE for iteration ' num2str(icall) ' is ' num2str(f)])
elseif strcmp(error_metric, 'RMSE')
  f = myRMSE(Qobs, Q);
  disp(['RMSE for iteration ' num2str(icall) ' is ' num2str(f)])
end

return
