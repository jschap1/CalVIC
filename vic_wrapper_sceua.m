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

% These parameters might be better placed in the cv_params file.
soil_line = 154;
log_out_line = 182;
out_line = 183;

icall = icall + 1; % using 1-indexing instead of 0-indexing

if control_params.all_outputs
    % Creates a new directory for each set of outputs
    control_params.vic_out_dir = [control_params.vic_out_dir num2str(icall)];
    control_params.glacier_out_dir = [control_params.glacier_out_dir num2str(icall)];
    control_params.rout_in_dir = [control_params.rout_in_dir num2str(icall)];
    control_params.rout_out_dir = [control_params.rout_out_dir num2str(icall)];
else
    % Replaces the contents of the previous directory with the new ones
    % Implementation is a bit clunky; could probably be done more
    % efficiently
    if icall > 1        
        disp('Implementation pending')
    end
    if icall==1
        control_params.vic_out_dir = [control_params.vic_out_dir num2str(1)];
        control_params.glacier_out_dir = [control_params.glacier_out_dir num2str(1)];
        control_params.rout_in_dir = [control_params.rout_in_dir num2str(1)];
        control_params.rout_out_dir = [control_params.rout_out_dir num2str(1)];
    end
end

% Create directory for VIC model outputs

mkdir(control_params.vic_out_dir)
disp(['Created directory ', control_params.vic_out_dir])

% Create directory for glacier melt model outputs
if control_params.add_glaciers == 1
    mkdir(control_params.glacier_out_dir)
    disp(['Created directory ', control_params.glacier_out_dir])
end

% Create directory for routing model inputs
mkdir(control_params.rout_in_dir)
disp(['Created directory ', control_params.rout_in_dir])

% Create directory for routing model outputs
mkdir(control_params.rout_out_dir)
disp(['Created directory ', control_params.rout_out_dir])

% Create directory for meta outputs
if icall==1
    mkdir(control_params.meta_output)
    disp(['Created directory ', control_params.meta_output])
end

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

n_pars = length(x); % number of calibration parameters

% Load soil parameters
if icall==1
     soils = load(control_params.init_soil_pars);
%      control_params.soil_param = control_params.init_soil_pars;
else
    soils = load(control_params.soil_param);
  
    % Replace the old parameters with the new parameters
    control_params.soil_pars % indices of soil parameters to calibrate
    for ii=1:n_pars
        soils(:, control_params.soil_pars(ii)) = x(ii);
    end

    % Write out the updated soil parameter file
    soil_path = fileparts(control_params.soil_param);
    control_params.soil_param = [soil_path '/soils_' num2str(icall) '.txt'];
    
    soil_format = '3l'; % 3l
    write_soils(5, soils, control_params.soil_param, soil_format);
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
A{soil_line} = ['SOIL ' control_params.soil_param];

% Change output directory
A{log_out_line} = ['LOG_DIR ' control_params.vic_out_dir '/'];
A{out_line} = ['RESULT_DIR ' control_params.vic_out_dir '/'];

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
    if control_params.use_hoffman == 1

        vic_exec_file = set_up_parallel(control_params);
        vic_run_command = ['qsub -cwd -V -N ' control_params.jobname ' -l h_data=' num2str(control_params.data_per_job) ...
            'M,h_rt=' num2str(control_params.time_per_job) ' -m n -t 1-' num2str(control_params.n_proc) ':1 ' vic_exec_file];
        system(vic_run_command);
        
        % Wait for VIC to finish executing before continuing
        while length(dir([control_params.vic_out_dir '/fluxes_*'])) < control_params.ncells
            pause(10)
        end

        disp('All files have been generated')
        disp(['Waiting ' num2str(control_params.wait_time) ' minutes to ensure all outputs are done'])
        pause(control_params.wait_time*60) % enter wait_time in minutes    

        % Just in case, check that all files have been generated and completed
        fnames = dir([control_params.vic_out_dir '/fluxes_*']);
        length(fnames)
        file_sizes = zeros(control_params.ncells, 1);
        for ii=1:control_params.ncells
            file_sizes(ii) = fnames(ii).bytes;
        end
        if min(file_sizes) == 0
            disp('Min file size is zero. Waiting for more time to ensure VIC model runs to completion.')
            disp(['Waiting ' num2str(control_params.wait_time) ' minutes to ensure all outputs are done'])
            pause(control_params.wait_time*60)          
        end        
        
    else    
    
        vic_run_command = [control_params.vic_command ' -g ' control_params.global_param_file];
        system(vic_run_command);
            
    end
    
end

% Check that all the files are done (not just started)
% fluxnames = dir([control_params.vic_out_dir '/fluxes_*']);
% nfiles = length(fluxnames);
% filesizes = zeros(nfiles, 1);
% for k=1:nfiles
%     thisfile = dir(fullfile(control_params.vic_out_dir, fluxnames(k).name));
%     filesizes(k) = thisfile.bytes;
% end

% This is a hacky solution that requires having a pretty good idea of how
% long it takes to process one grid cell. The amount of time to wait should
% be greater than the time to process one grid cell, just in case.
%
% The better solution would be to have VIC generate a special output
% (file?) when it finishes running, and just check for that, instead.

%% Adjust runoff for glacier melt contribution and write routing model input files

if control_params.lumped
%    disp('Adding glacier melt contribution');
%    disp('Glacier melt contribution is not currently supported for lumped modeling');
   add_glacier_contribution_lumped(control_params)
else
    % Check if outputs already exist (to make it easy to resume an interrupted run)
    fn1 = dir(fullfile(control_params.rout_in_dir, 'fluxes_*'));
    if isempty(fn1)
        ddf = control_params.ddf; % degree-day factor (mm/K/day)
        plotflag = 0;
        add_glacier_contribution(control_params.add_glaciers, ddf, ...
            control_params.glacier_fract_map, ...
            control_params.vic_out_dir, control_params.rout_in_dir, ...
            control_params.glacier_out_dir, plotflag)
    end
end

%% Routing model

if ~control_params.lumped
    
    % Write routing model input file
    B = read_global_param_file(control_params.rout_param_file);
    B{19} = fullfile(control_params.rout_in_dir, 'fluxes_');
    B{22} = [control_params.rout_out_dir '/'];
    write_global_param_file(B, control_params.rout_param_file)

    % Run routing model
    rout_run_command = [control_params.rout_command ' ' control_params.rout_param_file];
    disp('Running routing model')
    disp(rout_run_command)
    system(rout_run_command)

    % Remove existing unit hydrograph files
    !rm *.uh_s
    
end

%% Compute error metric

% Load gage discharge

monthly_flag = 1;
if monthly_flag
    disp('Monthly gage data')
    obs = readmatrix(control_params.discharge_obs); % default units are cms
    time_obs = datetime(obs(:,1), obs(:,2), 15);
else
    disp('Daily gage data')
    obs = dlmread(control_params.discharge_obs, '\t', 1, 0); % default units are cfs
    time_obs = datetime(obs(:,1), obs(:,2), obs(:,3));    
end
Qobs = obs(:,4);
Qobs(Qobs<0) = NaN; % remove negative "observed" values of discharge
% figure, plot(time_obs, Qobs)

if control_params.lumped
    % Workflow for lumped basin modeling/comparison
%     Qobs = (12/39.37)^3*Qobs; % convert to m3/s
%     disp('Converted discharge observations from cfs to m3s')
    graph_units = 'm^3/s';
    
    if control_params.add_glaciers
        tempname = dir(fullfile(control_params.rout_in_dir, 'fluxes_*'));
        lumped_outname = fullfile(control_params.rout_in_dir, tempname(1).name);
        rout = dlmread(lumped_outname, '\t', 3, 0);
        Q = rout(:,6) + rout(:,7);        
    else
        tempname = dir(fullfile(control_params.vic_out_dir, 'fluxes_*.txt'));
        lumped_outname = fullfile(control_params.vic_out_dir, tempname(1).name);
        rout = dlmread(lumped_outname, '\t', 3, 0);
        Q = rout(:,4) + rout(:,5);        
    end
    
%     rout = dlmread('/Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/out/vic_out_1/fluxes_31.91625_77.28219.txt', '\t', 3, 0);
    A = control_params.basin_area; % basin area, in sq km
    Q = Q*A*1000/(24*3600); 
    
% Examining the streamflow components
% time_rout = datetime(rout(:,1), rout(:,2), rout(:,3));
% Qb = A*rout(:,5)*1000/(24*3600); % baseflow (m3/s)
% Qd = A*rout(:,4)*1000/(24*3600); % runoff (m3/s)
% figure, subplot(3,1,1)
% plot(time_rout, Qb, 'linewidth', 2)
% title('Baseflow')
% grid on
% ylim([0,max(Qd+Qb)])
% subplot(3,1,2)
% plot(time_rout, Qd, 'linewidth', 2)
% title('Runoff')
% grid on
% ylim([0,max(Qd+Qb)])
% subplot(3,1,3)
% plot(time_rout, Qd+Qb, 'linewidth', 2)
% hold on
% title('Total Streamflow')
% grid on
% plot(time_obs, Qobs, 'linewidth', 2)
% legend('Sim','Obs')
% ylim([0,max(Qd+Qb)])
  
else
    graph_units = 'cfs';
    rout = load(fullfile(control_params.rout_out_dir, 'STA1 .day'));
    Q = rout(:,4);    
end

time_rout = datetime(rout(:,1), rout(:,2), rout(:,3));

% Note that the above line will fail if there are "NA" values in the
% observations. Instead, use negative values as the NAflag.

% obs = dlmread(control_params.discharge_obs, '\t', 3, 0);

% obs = dlmread(filenames.obs_discharge, '\t', 1, 0);
% Qobs = obs(:,4)/1000; % cfs (dividing by 1000 for the one-cell test)
% Qobs = obs(:,4) + obs(:,5); % combining runoff and baseflow for the one-cell test
% Qobs = Qobs + 10; % add a bias, just because

% Load discharge calculated with routing model
% rout = dlmread(['/Volumes/HD3/SWOTDA/Calibration/Pandoh/vic/out_' num2str(icall) '/fluxes_31.84375_77.03125.txt'], '\t', 3, 0);
% rout = dlmread(['/Volumes/HD3/SWOTDA/Calibration/onecell/vic_out_' num2str(icall) '/fluxes_31.84375_77.03125.txt'], '\t', 3, 0);
% Q = rout(:,4) + rout(:,5);

% Make sure the two time series are consistent
if monthly_flag
    [time_rout, Q] = daily_to_monthly(time_rout, Q, 'mean');
end

start_time = max(time_rout(1), time_obs(1));
end_time = min(time_rout(end), time_obs(end));
Qobs = Qobs(find(time_obs == start_time):find(time_obs == end_time));
Q = Q(find(time_rout == start_time):find(time_rout == end_time));
time_merged = time_rout(find(time_rout == start_time):find(time_rout == end_time));

% Plot once every 10 function calls
% if mod(icall,10) == 0
%     plotflag = 1;
% else
%     plotflag = 0;
% end
plotflag = 1;
if plotflag
    figure('Visible','off')
    set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')')
    plot(time_merged, Q)
    xlabel('Time')
    ylabel('Flow')
    title(['Discharge (' graph_units ')'])
    set(gca, 'fontsize', 18)   
    hold on
    plot(time_merged, Qobs)
    legend('Predicted','Observed','Location','NW')
    saveas(gcf, fullfile(control_params.rout_out_dir, ['discharge_plot_iter_' num2str(icall), '.png']))
    close(gcf)
end

if length(Q) ~= length(Qobs)
  error('Lengths of discharge records are not equal')
end

%% Calculate error metric, ignore the first few time steps (spin-up)

Qobs = Qobs(control_params.n_spinup:end);
Q = Q(control_params.n_spinup:end);
time_merged = time_merged(control_params.n_spinup:end);

error_metric = control_params.objective;
if strcmp(error_metric, 'NSE')
  f = -1*myNSE(Qobs, Q); % to maximize NSE, minimize -NSE
  f_alt = myRMSE(Qobs, Q); % calculate both parameters, either way
  disp(['NSE for iteration ' num2str(icall) ' is ' num2str(-1*f)])
elseif strcmp(error_metric, 'RMSE')
  f = myRMSE(Qobs, Q);
  f_alt = myNSE(Qobs, Q);
  disp(['RMSE for iteration ' num2str(icall) ' is ' num2str(f)])
elseif strcmp(error_metric, 'KGE')
  f = -1*myKGE(Qobs, Q);
  f_alt = myNSE(Qobs, Q);
  disp(['KGE for iteration ' num2str(icall) ' is ' num2str(-1*f)])
end

meta_outname = fullfile(control_params.meta_output, ['meta_output_' num2str(icall), '.mat']);
save(meta_outname, 'f','f_alt','Q','Qobs','time_merged','x')
% disp(['Saved metadata to ' meta_outname])

% pause

return
