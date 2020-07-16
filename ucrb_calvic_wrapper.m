% Updated version of CalVIC for the UCRB model
%
% July 16, 2020 JRS
%
%

addpath(genpath('/home/jschap/Documents/Codes/CalVIC'))
addpath(genpath('/home/jschap/Documents/Codes/VICMATLAB/vicmatlab'))

%% SCE-UA hyperparameters

maxn = 25; % maximum number of iterations
kstop = 15; % kstop (number of shuffling loops in which the criterion value must change by the given percentage before optimization is terminated)
pcento = 1; % pcento (percentage by which the criterion value must change in given number (kstop) of shuffling loops to continue optimization)
peps = 0.01; % convergence level for parameter set (lower number means smaller difference between parameters of the population required for stop)
iseed = 704753262; % initial (random) seed
iniflg = 1; % flag for whether to count initial iteration
ngs = 2; % number of complexes

%% Initial guess for parameters

% b, ds, dsmax, ws, t3, bd
% x0 = [0.26, 0.00392, 18.07, 0.575, 1.855, 1079];
% % x0 = [0.086, 0.34, 14.8, 0.748, 2.09, 1064];
% bl = [0.05, 0.01, 0.01, 0.5, 0.5, 1000]; % lower bound
% bu = [0.4, 0.5, 30, 1, 2.5, 1800]; % upper bound

% b, t3, dsmax
% x0 = [0.2438722, 2.094797, 28.66799];
x0 = [0.11, 1.8, 1.8];
% x0 = [0.1366702, 2.0151629, 1.6731925];
bl = [0.001, 0.5, 0.001];
bu = [0.4, 2.5, 30];

%% Call SCE-UA

% xf is the vector of objective function values
[bestx,bestf,xf] = sceua(x0,bl,bu,maxn,kstop,pcento,peps,ngs,iseed,iniflg);

%% Run one more time with optimal parameters

% [kge, t_m, q_m, t, q] = vic_wrapper_sceua([], x0, 9999, 1);
[kge, t_m, q_m, t, q] = vic_wrapper_sceua([], bestx, 9999, 1);

b_1 = 0.2438722;
t3_1 = 2.094797;
dsmax_1 = 28.66799;
x1 = [b_1,t3_1,dsmax_1];

% x1 = [0.00118, 0.00139, 7.67, 0.59, 1.77, 1201]; % run with dds parameters
% x1 = [0.2, 0.001, 19.6419, 0.9, 0.7, 1035.33]; % run with dds parameters
% x1 = [0.26, 0.0039, 18.07, 0.575, 1.86, 1079];
[kge, t_m, q_m, t, q] = vic_wrapper_sceua([], x1, 10001, 1);

writetable(table(q_m, q), './qtable.txt')
myKGE(q, q_m)

%% Make some plots for interpretation

figure, plot(t_m, q_m)
hold on, plot(t,q)
legend('pred','obs')
set(gca, 'fontsize', 16)

figure, plot(fliplr(-xf), 'o')
xlabel('iteration'), ylabel('kge')
grid on
set(gca, 'fontsize', 16)

figure
plot(q,q_m, '.')
loglog(q, q_m, '.')
xlabel('observed'), ylabel('predicted')
grid on
refline(1,0)
axis equal
axis tight
set(gca, 'fontsize', 16)
