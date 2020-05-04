% Adds glacier contribution to runoff
%
% Uses degree day modeling approach
% Use this to update the VIC flux files with the adjusted runoff
%
% INPUTS
% ddf = degree day factor for snow, mm/deg. C/day
% glacier_fraction = glacier fraction raster
% indir = directory where flux files are stored
% outdir = directory to write input files for the routing model
% saveloc = directory to save data as a .mat file
% plotflag = flag for plotting
%
% TODO
% It might be appropriate to also remove the portion of baseflow that
% originates in the glaciated fraction of a pixel

function add_glacier_contribution_lumped(control_params)

%% Initialization

disp('Modifying runoff to account for glacier contribution')
disp(['Degree-day factor assumed to be ' num2str(control_params.ddf) ' mm/K/day'])

% load surface temperature and runoff data
fluxnames = dir(fullfile(control_params.vic_out_dir, 'fluxes*'));

raw_output = dlmread(fullfile(control_params.vic_out_dir, fluxnames(1).name), '\t', 3, 0);
nsteps = size(raw_output, 1);
ncells = length(fluxnames);

disp(['Modifying runoff for ' num2str(ncells) ' grid cells'])

runoff_all = zeros(ncells, nsteps);
baseflow_all = zeros(ncells, nsteps);
temperature_all = zeros(ncells, nsteps);
lat_all = zeros(ncells, 1);
lon_all = zeros(ncells, 1);

disp('Loading runoff, baseflow, and temperature data from VIC results')
for k=1:ncells

    % get lat/lon for the cell
    tmpstring = fluxnames(k).name;
    tmpstring = strrep(tmpstring,'-','');
    tmpstring = strrep(tmpstring,'.','_');
    tmpstring = strrep(tmpstring,'fluxes_','');
    tmpstring = strrep(tmpstring,'_txt','');
    tmpcell = strsplit(tmpstring,'_');
    lat_all(k) = str2double([tmpcell{1} '.' tmpcell{2}]);
    lon_all(k) = str2double([tmpcell{3} '.' tmpcell{4}]);

    raw_output = dlmread(fullfile(control_params.vic_out_dir, fluxnames(k).name), '\t', 3, 0);

    runoff_all(k,:) = raw_output(:,4);
    baseflow_all(k,:) = raw_output(:,5);
    temperature_all(k,:) = raw_output(:,6);
    
    if mod(k, 1e3) == 0
        disp(['Progress: ' num2str(k,2) ' out of ' num2str(ncells)])
    end

end

%% Recalculate total runoff

if control_params.add_glaciers
    
    % load glacier fraction data
    glacierfract = load(control_params.glacier_fract_map);
    gf = glacierfract.glacier_fract;
    % Should use a different glacier fraction for each timestep

    %%% To turn off the glacier contribution ----------------------------
    % glacierfract = zeros(size(glacierfract));
    % -------------------------------------------------------------------
    
    % Calculate melt runoff using degree-day model
    PDD = temperature_all;
    PDD(temperature_all<0) = 0;
    m = control_params.ddf*PDD; % melt flux (mm)    
    
    % Note that this syntax is valid beginning with Matlab R2016, when
    % Mathworks introduced arithmetic expansion. For earlier versions, a
    % different formula is required.
    % r_total = gf.*m + (1-gf).*runoff_all;
    sz2 = size(m, 2);
    r_total = repmat(gf, 1, sz2).*m + (1-repmat(gf, 1, sz2)).*runoff_all;
    
    disp('Removing baseflow from glaciated area')
    disp('This may/may not be the correct decision')
    baseflow_total = (1-repmat(gf, 1, sz2)).*baseflow_all;
    % it doesn't seem to make a big difference, at least for Pandoh test
    % case, anyway
    
    % test
    % gf = rand(134,1);
    % m = 10*ones(134,31);
    % runoff_all = 15*ones(134,31);
    % r_total_2 = repmat(gf, 1, 31).*m + (1-repmat(gf, 1, 31)).*runoff_all;
    % ma_expansion = repmat(ma,3,1)

    save(fullfile(control_params.glacier_out_dir, 'glacier_contribution.mat'), 'r_total', 'gf', 'runoff_all', 'lat_all', 'lon_all')
    disp(['Saved glacier melt contribution to ' control_params.glacier_out_dir])

end % accounting for glaciers

%% Write modified runoff to routing model input files

disp(['Writing routing model input files to ' control_params.rout_in_dir]);

for k=1:ncells
    if control_params.add_glaciers
        routing_input = [raw_output(:,1:3) zeros(nsteps, 2) r_total(k,:)' baseflow_total(k,:)'];
    else
        routing_input = [raw_output(:,1:3) zeros(nsteps, 2) runoff_all(k,:)' baseflow_all(k,:)'];
    end

    dlmwrite(fullfile(control_params.rout_in_dir, fluxnames(k).name(1:end-4)), routing_input, '\t'); % removes .txt suffix
end
