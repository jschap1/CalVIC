% Post-processing
%
% Make plots from calibration output

clear

%%
metadir = '/Volumes/HD3/SWOTDA/Calibration/Pandoh/out/meta';
figdir = '/Volumes/HD3/SWOTDA/Calibration/Pandoh/Figures';
mkdir(figdir)
fnames = dir([metadir '/*.mat']);
ncalls = length(fnames); % number of function calls

%%

rmse = zeros(ncalls, 1);
nse =  zeros(ncalls, 1);
for k=1:ncalls
    load(fullfile(metadir, fnames(k).name))
    rmse(k) = f;
    nse(k) = f_alt;
end

figure, subplot(2,1,1), jsplot(1:ncalls, rmse, '', 'Function call', 'RMSE (cfs)', 18)
subplot(2,1,2), jsplot(1:ncalls, nse, '', 'Function call', 'NSE (-)', 18)

%%
% First iteration
icall = 1;
load(fullfile(metadir, fnames(icall).name))
figure
jsplot(time_merged, Q, 'Flow', 'Time', 'Discharge (cfs)', 18)
hold on
plot(time_merged, Qobs)
legend('Predicted','Observed','Location','NW')
saveas(gcf, fullfile(figdir, ['discharge_plot_iter_' num2str(icall), '.png']))
    
%% Best iteration

[rmse_min, min_ind] = min(rmse);
icall = min_ind;
load(fullfile(metadir, fnames(icall).name))
figure
jsplot(time_merged, Q, 'Flow', 'Time', 'Discharge (cfs)', 18)
hold on
plot(time_merged, Qobs)
legend('Predicted','Observed','Location','NW')
% saveas(gcf, fullfile(figdir, ['discharge_plot_iter_' num2str(icall), '.png']))

%% Both on same axis

figure
load(fullfile(metadir, fnames(1).name))
jsplot(time_merged, Q, 'Flow', 'Time', 'Discharge (cfs)', 18)
hold on
plot(time_merged, Qobs)
load(fullfile(metadir, fnames(min_ind).name))
plot(time_merged, Q)
legend('Predicted (First)','Observed', 'Predicted (Best)', 'Location','NW')

%% Plot for IDRE poster

% figdir = '/Volumes/HD3/SWOTDA/Presentations/IDRE_Figures';

% ff = 35.314666212661;
ff = 1; % optional conversion factor

% First iteration
icall = 1;
load(fullfile(metadir, fnames(icall).name))
figure
plot(time_merged, Qobs/ff, 'linewidth', 3, 'color', 'black')
xlabel('Time')
ylabel('Discharge (cfs)')
title('Pandoh Dam')
hold on
plot(time_merged, Q/ff, 'linewidth', 3, 'color', 'blue')
load(fullfile(metadir, fnames(min_ind).name))
plot(time_merged, Q/ff, 'linewidth', 3, 'color', 'green')
legend('Observed','Predicted (initial)', 'Predicted (calibrated)', 'Location','NW')
set(gca, 'fontsize', 36)

%% Sensitivity analysis

metadir = '/Volumes/HD3/SWOTDA/Calibration/Pandoh/out/meta';
% metadir = '/Volumes/HD3/SWOTDA/Calibration/TuoSub/out/meta';
fnames = dir([metadir '/*.mat']);
n = length(fnames);
A = cell(n,1);
for i=1:n
    A{i} = load(fullfile(metadir, fnames(i).name));
end

% Plot the response of RMSE to each variable
rmse = zeros(n,1);
b = zeros(n,1);
Ds = zeros(n,1);
Dsmax = zeros(n,1);
Ws = zeros(n,1);
Qobs = A{1}.Qobs;
time_merged = A{1}.time_merged;
nt = length(Qobs);
Q = zeros(nt,n);
for i=1:n
    rmse(i) = A{i}.f;
    b(i) = A{i}.x(1);
    Ds(i) = A{i}.x(2);
    Dsmax(i) = A{i}.x(3);
    Ws(i) = A{i}.x(4);
    Q(:,i) = A{i}.Q;
end

[min_rmse, min_ind] = min(rmse);

%%
figure
msize = 25;

subplot(2,2,1)
scatter(b, rmse, msize, 'filled')
hold on
plot(b(min_ind), min_rmse, 'r.', 'markersize', msize)
xlabel('b (-)')
ylabel('RMSE (cfs)')
grid on
set(gca, 'fontsize', 16)

subplot(2,2,2)
scatter(Ds, rmse, msize, 'filled')
hold on
plot(Ds(min_ind), min_rmse, 'r.', 'markersize', msize)
xlabel('Ds (-)')
ylabel('RMSE (cfs)')
grid on
set(gca, 'fontsize', 16)

subplot(2,2,3)
scatter(Dsmax, rmse, msize, 'filled')
hold on
plot(Dsmax(min_ind), min_rmse, 'r.', 'markersize', msize)
xlim([0,30])
xlabel('Dsmax (mm/day)')
ylabel('RMSE (cfs)')
grid on
set(gca, 'fontsize', 16)

subplot(2,2,4)
scatter(Ws, rmse, msize, 'filled')
hold on
plot(Ws(min_ind), min_rmse, 'r.', 'markersize', msize)
xlabel('Ws (-)')
ylabel('RMSE (cfs)')
grid on
set(gca, 'fontsize', 16)


%% Plot all discharge series

figure, 
hold on
plot(time_merged, Qobs, 'black', 'linewidth', 3)
xlabel('T')
ylabel('Q (cfs)')
for i=1:n
    plot(time_merged, Q(:,i), 'cyan', 'linewidth', 1)
end
plot(time_merged, Q(:,min_ind), 'red', 'linewidth', 2)
plot(time_merged, Q(:,1), 'blue', 'linewidth', 1)
set(gca, 'fontsize', 16)


