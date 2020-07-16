% VIC wrapper for the shuffled complex algorithm
% October 31, 2019 JRS
%
% Updated to customize for UCRB calibration 7/6/2020 JRS
%
% INPUTS
% x = parameters
% icall = number of the current model run
%
% OUTPUTS
% f = error metric (NSE)

function [f, t_m, q_m, t, q] = vic_wrapper_sceua(nopt, x, icall, plotflag)

icall = icall + 1; % using 1-indexing instead of 0-indexing

%% Write soil parameter file using the parameters x

spf = '/home/jschap/Documents/ESSD/lumped_cal/sceua/soils.txt';
soils = load(spf);

soils(:,5) = x(:,1); % b
soils(:,25) = x(:,2); % t3
soils(:,7) = x(:,3); % dsmax

% soils(:,8) = x(:,4); % ws
% soils(:,6) = x(:,2); % ds
% soils(:,34) = x(:,6); % bd
% soils(:,35) = x(:,6); % bd
% soils(:,36) = x(:,6); % bd

write_soils(5, soils, spf, '3l');

%% Run VIC

gpf = '/home/jschap/Documents/ESSD/lumped_cal/sceua/global_param.txt';
vic = '/home/jschap/Documents/Software/VIC/vic/drivers/classic/vic_classic.exe';
vic_run_command = [vic ' -g ' gpf];
system(vic_run_command);            

%% Calculate summed runoff and baseflow

wbfile = '/home/jschap/Documents/ESSD/lumped_cal/dds/WY2002-2005_kitchen_sink/validation/wb_39.24321_-109.02250.txt';
wb = readmatrix(wbfile);
nspinup = 365;
wb = wb((nspinup+1):end, :);

qd = wb(:,5);
qb = wb(:,6);
A_basin = 293600; % km2, size of basin
conversion_factor = 1000*A_basin/(24*3600);
q_hat = conversion_factor*(qb+qd);
t_hat = datetime(wb(:,1), wb(:,2), wb(:,3));
% figure, plot(t_hat, q_hat)

% Calculate monthly flow
[t_m, q_m] = daily_to_monthly(t_hat, q_hat, 'mean');
% figure, plot(t_m, q_m)

%% Load in validation data
  
valdata = '/hdd/ESSD/data/naturalized_flow_usgs_09380000.txt';
truth = readmatrix(valdata);
qobs = truth(:,4);
tobs = datetime(truth(:,1), truth(:,2), 15);

[~, i1] = ismember(t_m(1), tobs);
[~, i2] = ismember(t_m(end), tobs);

t = tobs(i1:i2);
q = qobs(i1:i2);

if plotflag
    lw = 1.5;
    figure, plot(t_m, q_m, 'linewidth', lw)
    hold on, plot(t,q, 'linewidth', lw)
    xline(datetime(2005, 10, 1), '--', 'linewidth', lw)
    legend('pred','obs', 'cal/val')   
    xlabel('Time')
    ylabel('Monthly streamflow (m^3/s)')
    set(gca, 'fontsize', 16)
end

%% Compute error metric

error_metric = 'KGE';
if strcmp(error_metric, 'NSE')
  f = -1*myNSE(q, q_m); % to maximize NSE, minimize -NSE
  disp(['NSE for iteration ' num2str(icall) ' is ' num2str(-1*f)])
elseif strcmp(error_metric, 'RMSE')
  f = myRMSE(q, q_m);
  disp(['RMSE for iteration ' num2str(icall) ' is ' num2str(f)])
elseif strcmp(error_metric, 'KGE')
  f = -1*myKGE(q, q_m);
  disp(['KGE for iteration ' num2str(icall) ' is ' num2str(-1*f)])
end


return
