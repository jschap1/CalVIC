% Post-processing
%
% Make plots from calibration output

clear

metadir = '/Volumes/HD3/SWOTDA/Calibration/Pandoh2/meta';
figdir = '/Volumes/HD3/SWOTDA/Calibration/Pandoh2/Figures';
mkdir(figdir)
fnames = dir([metadir '/*.mat']);
ncalls = length(fnames); % number of function calls

rmse = zeros(ncalls, 1);
nse =  zeros(ncalls, 1);
for k=1:ncalls
    load(fullfile(metadir, fnames(k).name))
    rmse(k) = f;
    nse(k) = f_alt;
end

figure, jsplot(1:ncalls, rmse, '', 'Function call', 'RMSE (cfs)', 18)
figure, jsplot(1:ncalls, nse, '', 'Function call', 'NSE (-)', 18)

% First iteration
icall = 1;
load(fullfile(metadir, fnames(icall).name))
figure
jsplot(time_merged, Q, 'Flow', 'Time', 'Discharge (cfs)', 18)
hold on
plot(time_merged, Qobs)
legend('Predicted','Observed','Location','NW')
saveas(gcf, fullfile(figdir, ['discharge_plot_iter_' num2str(icall), '.png']))
    
%% Last iteration
icall = 34;
load(fullfile(metadir, fnames(icall).name))
figure
jsplot(time_merged, Q, 'Flow', 'Time', 'Discharge (cfs)', 18)
hold on
plot(time_merged, Qobs)
legend('Predicted','Observed','Location','NW')
saveas(gcf, fullfile(figdir, ['discharge_plot_iter_' num2str(icall), '.png']))

%% Both on same axis

figure
load(fullfile(metadir, fnames(1).name))
jsplot(time_merged, Q, 'Flow', 'Time', 'Discharge (cfs)', 18)
hold on
plot(time_merged, Qobs)
load(fullfile(metadir, fnames(34).name))
plot(time_merged, Q)
legend('Predicted (1)','Observed', 'Predicted (34)', 'Location','NW')


