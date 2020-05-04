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
% parameter_file = '/Volumes/HD3/SWOTDA/Calibration/Lumped/Mangla/cv_params.sh';
% parameter_file = '/Volumes/HD3/SWOTDA/Calibration/Lumped/Tarbela/cv_params.sh';
% parameter_file = '/Volumes/HD3/SWOTDA/Calibration/Lumped/Bhakra/cv_params.sh';
% parameter_file = '/Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/cv_params.sh';

B = read_parfile(parameter_file);

params = struct();

% Locations of software packages
params.vic_command = strsplit_spec(B{1});
params.rout_command = strsplit_spec(B{2});
params.sce_dir = strsplit_spec(B{3});

% Input files
params.global_param_file = strsplit_spec(B{4});
params.init_soil_pars = strsplit_spec(B{5});
params.rout_param_file = strsplit_spec(B{6}); % note that the stnloc station must be named STA1
params.discharge_obs = strsplit_spec(B{7});
params.glacier_fract_map = strsplit_spec(B{8});

% Control parameters
params.vic_out_dir = strsplit_spec(B{9});
params.glacier_out_dir = strsplit_spec(B{10});
params.rout_in_dir = strsplit_spec(B{11});
params.rout_out_dir = strsplit_spec(B{12});

% SCE-UA parameters
params.n_spinup = str2double(strsplit_spec(B{13}));
params.n_complexes = str2double(strsplit_spec(B{14}));
params.max_iter = str2double(strsplit_spec(B{15})); % maximum allowable number of function evaluations
params.kstop = str2double(strsplit_spec(B{16})); % maximum number of evolution loops before convergence
params.pcento = str2double(strsplit_spec(B{17})); % percentage change allowed in kstop loops before convergence
params.peps = str2double(strsplit_spec(B{18})); % controls when the SCE method converges
params.iseed = str2double(strsplit_spec(B{19})); % random seed number for generating initial guesses
params.iniflg = str2double(strsplit_spec(B{20})); % flag for initial parameter array
params.soil_pars = str2double(strsplit(strsplit_spec(B{21}), ',')); % parameters to calibrate
params.par_bounds = strsplit_spec(B{22}); % upper and lower bounds for calibration parameters

% Parallelization parameters
params.use_hoffman = str2double(strsplit_spec(B{23}));
params.jobname = strsplit_spec(B{24});
params.data_per_job = str2double(strsplit_spec(B{25})); % MB
params.time_per_job = str2double(strsplit_spec(B{26})); % hours
params.n_proc = str2double(strsplit_spec(B{27})); % number of processors to use
params.wait_time = str2double(strsplit_spec(B{28})); % minutes
params.meta_output = strsplit_spec(B{29}); % directory for meta outputs, like RMSE(t), etc.
% time to wait in between runs 
% (ensures that VIC files are completely done processing before moving on with the calibration

% Glacier parameters
params.add_glaciers = str2double(strsplit_spec(B{30})); % flag (0 or 1)
params.ddf = str2double(strsplit_spec(B{31})); % mm/K/day;

% Other parameters
params.objective = strsplit_spec(B{32}); % RMSE or NSE
params.all_outputs = str2double(strsplit_spec(B{33}));
params.lumped = str2double(strsplit_spec(B{34})); % Lumped or distributed
params.basin_area = str2double(strsplit_spec(B{35})); % km2, used for lumped model
% addpath(params.sce_dir);
% addpath('/Users/jschap/Documents/MATLAB/sce_matlab');

%%
% ------------------------------------------------------------------------

% clearvars -except soils_vg;
global BESTX BESTF ICALL PX PF

ngs = params.n_complexes; % number of complexes
% ngs = 6; % more complexes, better chance of finding a global optimum
% ngs = 2; % few complexes, fewer function evaluations may be required

% soils_vg = load('/Volumes/HD3/VICParametersGlobal/Global_1_16/v1_4/Classic/soils_3L_MERIT.txt');
% extent = '/Volumes/HD3/SWOTDA/FDT/v10282019/pandoh_basinmask_coarse.tif';
% extent = '/Volumes/HD3/SWOTDA/FDT/v10282019/tarbela_basinmask_coarse.tif';
% soils_sub = subset_soils(soils_vg, extent, '/Volumes/HD3/SWOTDA/Calibration/Tarbela/soils_tarbela.txt', '3l', 5);
% soils_pandoh = load('/Volumes/HD3/SWOTDA/Calibration/Pandoh/soils_pandoh.txt');

% initial soil parameter file
soils_init = load(params.init_soil_pars);

ncells = size(soils_init, 1); % should be ~134 for Pandoh Dam watershed
disp(['The number of grid cells is ' num2str(ncells)])

params.ncells = ncells; % for input into vic_wrapper

% Read in calibration parameter bounds
[bl, bu] = read_bounds(params.soil_pars, params.par_bounds);

% Get initial guesses
x0 = soils_init(1,params.soil_pars);

%% SCE-UA hyperparameters

maxn = params.max_iter;
kstop = params.kstop;
pcento = params.pcento;
peps = params.peps;
iseed = params.iseed;
iniflg = params.iniflg;

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

