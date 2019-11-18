% SCE-UA wrapper
%
% Runs the SCE-UA calibration algorithm
%
% Update 11/15/2019 JRS
% Rewritten as a function that can be compiled with mcc
% Designed to use a control parameter file to define file locations, etc.
%
% Update 11/18/2019 JRS
% Now written to make use of multiple processors by splitting the soil
% parameter file into multiple files

function exitcode = run_sceua(parameter_file)

% addpath('/Users/jschap/Documents/Codes/VICMATLAB/Subsetting')
% addpath('/Users/jschap/Documents/Codes/VICMATLAB/Make_Soils')
% addpath('/Users/jschap/Documents/Codes/VICMATLAB/Control')
% cd('/Volumes/HD3/SWOTDA/Calibration')
% addpath(genpath('/Volumes/HD3/SWOTDA/Calibration/CalVIC'))

%% Read the parameter file ------------------------------------------------

% parameter_file = '/Volumes/HD3/SWOTDA/Calibration/Pandoh2/cv_params.txt';

B = read_parfile(parameter_file);

params = struct();
params.vic_command = strsplit_spec(B{1});
params.rout_command = strsplit_spec(B{2});
params.sce_dir = strsplit_spec(B{3});

params.global_param_file = strsplit_spec(B{4});
params.init_soil_pars = strsplit_spec(B{5});
params.rout_param_file = strsplit_spec(B{6});
params.discharge_obs = strsplit_spec(B{7});
params.glacier_fract_map = strsplit_spec(B{8});

params.vic_out_dir = strsplit_spec(B{9});
params.glacier_out_dir = strsplit_spec(B{10});
params.rout_in_dir = strsplit_spec(B{11});
params.rout_out_dir = strsplit_spec(B{12});

params.n_spinup = str2double(strsplit_spec(B{13}));
params.n_complexes = str2double(strsplit_spec(B{14}));

params.n_proc = str2double(strsplit_spec(B{15})); % number of processors to use
params.jobarray_template = strsplit_spec(B{16}); 

addpath(params.sce_dir);
% addpath('/Users/jschap/Documents/MATLAB/sce_matlab');

% ------------------------------------------------------------------------

% clearvars -except soils_vg;
global BESTX BESTF ICALL PX PF

ngs = params.n_complexes; % number of complexes
% ngs = 6; % more complexes, better chance of finding a global optimum
% ngs = 2; % few complexes, fewer function evaluations may be required

% soils_vg = load('/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_4/Classic/soils_3L_MERIT.txt');
% extent = '/Volumes/HD3/SWOTDA/FDT/v10282019/pandoh_basinmask_coarse.tif';
% soils_pandoh = subset_soils(soils_vg, extent, './Pandoh/soils_pandoh.txt', '3l', 5);
% soils_pandoh = load('/Volumes/HD3/SWOTDA/Calibration/Pandoh/soils_pandoh.txt');

% initial soil parameter file
soils_init = load(params.init_soil_pars);

ncells = size(soils_init, 1); % should be ~134 for Pandoh Dam watershed
disp(['The number of grid cells is ' num2str(ncells)])

%%
% Calibration parameters (all uniform): 
% b_infilt [10^-5 to 0.4] % higher values increase runoff (column 5)
% ds [0.001 to 1] % fraction of the Dsmax parameter at which non-linear
% base-flow occurs (column 6)
% dsmax (column 7)
% ws [1e-4 to 1] % fraction of maximum soil moisture where non-linear
% baseflow occurs; usually is greater than 0.5 (column 8)
% layer 2, 3 soil depth % thicker soils decrease baseflow (columns 24, 25)

% bl = [1e-5, 1e-4, 0.01, 1e-4, 0.1, 0.3];
% bu = [0.4, 1, 20, 1, 0.3, 1.5];
% x0 = [soils_init(1,5), soils_init(1,6), soils_init(1,7), soils_init(1,8), soils_init(1,24), soils_init(1,25)];

% bl = 0.3*ones(1,ncells);
% bu = 1.5*ones(1,ncells);
% x0 = soils_pandoh(:,25)';

% % Spatially distributed b_infilt
% bl = 1e-5*ones(1,ncells);
% bu = 0.9*ones(1,ncells);
% x0 = soils_pandoh(:,5)';

% Everything but soil depth
% bl = [1e-5, 1e-4, 0.01, 1e-4];
% bu = [0.4, 1, 30, 1];
% x0 = [soils_init(1,5), soils_init(1,6), soils_init(1,7), soils_init(1,8)];

% Calibrating on spatially uniform b_infilt, Ds, and Ws, 
bl = [1e-5, 1e-4, 1e-4];
bu = [0.4, 1, 1];
x0 = [mean(soils_init(:,5)), mean(soils_init(:,6)), mean(soils_init(:,8))];

% Function for running the VIC model and computing an error metric
% !copy vic_wrapper_sceua.m functn.m
% Note: the ! means the command is to be executed by the OS
% It is shorthand for system()

%% SCE-UA hyperparameters

maxn = 50;
kstop = 10;
pcento = 1;
peps = 0.01;
iseed = -1;
iniflg = 1;

%% Do SCE-UA

[bestx,bestf] = sceua(params, x0,bl,bu,maxn,kstop,pcento,peps,ngs,iseed,iniflg);
 
% Keep in mind how long it takes to run one iteration of the VIC model
% For the Pandoh basin, with 134 grid cells, it takes about 4.5 hours to
% run VIC using energy balance, frozen soils, and snowbands for 2010-2019.
%
% In water balance mode only, it takes about 47 minutes for 2010-2019.
% In water balance mode only, it takes about 4 minutes for 2015.
% BUT --> it can get stuck when you change the parameters!

exitcode = 1;

%% Run one more simulation with the estimated parameters
% 
% icall_final = 9999;
% f = vic_wrapper_sceua(length(bestx), bestx, icall_final);

return

