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

function add_glacier_contribution(ddf, glacier_fraction, indir, outdir, saveloc, plotflag)

%% Initialization

disp('Modifying runoff to account for glacier contribution')
disp(['Degree-day factor assumed to be ' num2str(ddf) ' mm/K/day'])

if plotflag
  disp('Plots will be generated')
else
  disp('No plots will be generated')
end

% load glacier fraction data
[glacierfract, R] = geotiffread(glacier_fraction);
% Should use a different glacier fraction for each timestep

%%% To turn off the glacier contribution ----------------------------
% glacierfract = zeros(size(glacierfract));
% -------------------------------------------------------------------

% load surface temperature and runoff data
fluxnames = dir(fullfile(indir, 'fluxes*'));

raw_output = dlmread(fullfile(indir, fluxnames(1).name), '\t', 3, 0);
nsteps = size(raw_output, 1);
ncells = length(fluxnames);

disp(['Modifying runoff for ' num2str(ncells) ' grid cells'])

% these can be pretty big
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

    raw_output = dlmread(fullfile(indir, fluxnames(k).name), '\t', 3, 0);

    runoff_all(k,:) = raw_output(:,4);
    baseflow_all(k,:) = raw_output(:,5);
    temperature_all(k,:) = raw_output(:,6);
    
    if mod(k, 1e3) == 0
        disp(['Progress: ' num2str(k,2) ' out of ' num2str(ncells)])
    end

end

%% Calculate melt runoff using degree-day model

% get glacier fraction for this grid cell (location is inexact because the
% domains are slightly different)
[gc_r, gc_c] = latlon2pix(R, lat_all, lon_all);
gc_r = round(gc_r);
gc_c = round(gc_c);

gf = zeros(ncells, 1);
for k=1:ncells
    gf(k) = glacierfract(gc_r(k), gc_c(k));
end

PDD = temperature_all;
PDD(temperature_all<0) = 0;
m = ddf*PDD; % melt flux (mm)

%% Recalculate total runoff

% Note that this syntax is valid beginning with Matlab R2016, when
% Mathworks introduced arithmetic expansion. For earlier versions, a
% different formula is required.

% r_total = gf.*m + (1-gf).*runoff_all;
sz2 = size(m, 2);
r_total = repmat(gf, 1, sz2).*m + (1-repmat(gf, 1, sz2)).*runoff_all;

% test
% gf = rand(134,1);
% m = 10*ones(134,31);
% runoff_all = 15*ones(134,31);
% r_total_2 = repmat(gf, 1, 31).*m + (1-repmat(gf, 1, 31)).*runoff_all;
% ma_expansion = repmat(ma,3,1)

save(fullfile(saveloc, 'glacier_contribution.mat'), 'r_total', 'gf', 'runoff_all', 'lat_all', 'lon_all')
disp(['Saved glacier melt contribution to ' saveloc])

%% Write modified runoff to routing model input files

disp(['Writing routing model input files to ' outdir]);

account_for_glaciers = 1;
for k=1:ncells
    
    if account_for_glaciers
        routing_input = [raw_output(:,1:3) zeros(nsteps, 2) r_total(k,:)' baseflow_all(k,:)'];
    else
        routing_input = [raw_output(:,1:3) zeros(nsteps, 2) runoff_all(k,:)' baseflow_all(k,:)'];
    end

    dlmwrite(fullfile(outdir, fluxnames(k).name(1:end-4)), routing_input, '\t'); % removes .txt suffix
end

%% Make maps

if plotflag

  figsavedir = '/Volumes/HD3/SWOTDA/Calibration/Pandoh/Figures';

    [pandoh.dem, pandoh.R, pandoh.lon, pandoh.lat] = geotiffread2('/Volumes/HD3/SWOTDA/Calibration/Pandoh/pandoh_dem.tif');
    pandoh.dem(pandoh.dem==0) = NaN;
    figure
    plotraster(pandoh.lon, pandoh.lat, pandoh.dem, 'Pandoh Elevation (m)', 'Lon', 'Lat')
    pandoh.mask = pandoh.dem;
    pandoh.mask(~isnan(pandoh.mask)) = 1;
    saveas(gcf, fullfile(figsavedir, 'pandoh_dem.png'))

    glacierfract1 = flipud(geotiffread('/Volumes/HD3/SWOTDA/Calibration/Pandoh/pandoh_gf_crop.tif'));
    figure
    plotraster(pandoh.lon, pandoh.lat, glacierfract1.*pandoh.mask, 'Glacier fraction', 'Lon', 'Lat')
    saveas(gcf, fullfile(figsavedir, 'pandoh_glacier_fraction.png'))

    mean_glacier_contribution = mean(gf.*m,2);
    mean_glacier_contribution_map = xyz2grid(lon_all, lat_all, mean_glacier_contribution);
    figure
    plotraster(pandoh.lon, pandoh.lat, mean_glacier_contribution_map.*pandoh.mask, 'Glacier melt flux (mm)', 'Lon', 'Lat')
    saveas(gcf, fullfile(figsavedir, 'pandoh_glacier_melt.png'))

    mean_runoff_contribution = mean((1-gf).*runoff_all, 2);
    mean_runoff_contribution_map = xyz2grid(lon_all, lat_all, mean_runoff_contribution);
    figure
    plotraster(pandoh.lon, pandoh.lat, mean_runoff_contribution_map.*pandoh.mask, 'Runoff (mm)', 'Lon', 'Lat')
    saveas(gcf, fullfile(figsavedir, 'pandoh_runoff_og.png'))

    mean_temperature_contribution = mean(temperature_all, 2);
    mean_temperature_contribution_map = xyz2grid(lon_all, lat_all, mean_temperature_contribution);
    figure
    plotraster(pandoh.lon, pandoh.lat, mean_temperature_contribution_map.*pandoh.mask, 'Temperature (deg. C)', 'Lon', 'Lat')
    saveas(gcf, fullfile(figsavedir, 'pandoh_temperature_og.png'))

end

return

% %%
% % The old way (doesn't require as much RAM), but format is less useful
%
% for k=1:ncells
%
%     % get lat/lon for the cell
%     tmpstring = fluxnames(k).name;
%     tmpstring = strrep(tmpstring,'-','');
%     tmpstring = strrep(tmpstring,'.','_');
%     tmpstring = strrep(tmpstring,'fluxes_','');
%     tmpstring = strrep(tmpstring,'_txt','');
%     tmpcell = strsplit(tmpstring,'_');
%     lat = str2double([tmpcell{1} '.' tmpcell{2}]);
%     lon = str2double([tmpcell{3} '.' tmpcell{4}]);
%
%     raw_output = dlmread(fullfile(resultsdir, fluxnames(k).name), '\t', 3, 0);
%     temperature = raw_output(:,20);
%     runoff = raw_output(:,6);
%
%     % get glacier fraction for this grid cell (location is inexact because the
%     % domains are slightly different)
%     [gc_r, gc_c] = latlon2pix(R, lat, lon);
%     gc_r = round(gc_r);
%     gc_c = round(gc_c);
%     gf = glacierfract(gc_r, gc_c);
%
%     % calculate melt runoff using degree-day model
%     PDD = zeros(size(temperature));
%     PDD(temperature>0) = temperature(temperature>0); % positive degree-day sum (deg. C - day)
%     m = ddf*PDD; % melt flux (mm)
%
%     % Issue: temperatures are too cold. This may be due to the lapse rate.
%     % This is connected to the issue of insufficient snow melt-out
%
%     % recalculate total runoff
%     r_total = gf*m + (1-gf)*runoff;
%
%     % plot
%     figure, hold on
%     plot((1-gf)*runoff, 'green')
%     plot(gf*m, 'red')
%     plot(r_total, 'blue')
%     legend('No-glacier runoff','Melt runoff','Total runoff')
%
%     % Format for routing model input file:
%     % YYYY MM DD SKIP SKIP RUNOFF BASEFLOW
%     % raw_output(:,7) % baseflow
%     % raw_output(:,1:3) % year, month, day, skip, skip
%
%     % write modified runoff to routing model input files
%     routing_input = [raw_output(:,1:3) zeros(nsteps, 2) r_total raw_output(:,7)];
%     dlmwrite(fullfile(outdir, fluxnames(k).name), routing_input, '\t');
%
% end
