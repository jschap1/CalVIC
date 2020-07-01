% Distributed to lumped
%
% Makes a distributed VIC setup into a lumped VIC setup, so the model can
% be run quickly for calibration. Eventually, it would be cool to
% generalize this code to simply upscale the VIC setup from a finer
% resolution to a coarser resolution

clearvars -except soils

%% Inputs

basin = struct();
[basin.mask, basin.R, basin.lon, basin.lat] = geotiffread2('/home/jschap/Documents/ESSD/data/colo_mask.tif');
soilfile = '/home/jschap/Documents/ESSD/clipped_soils_VG.txt';
basin.soils = load(soilfile);
prefix = 'full_data_';
forcdir = '/media/jschap/HD_ExFAT/ucrb/disagg_forcings_L13';

% Lumped forcing outname
centroid = [-109.0225, 39.24321]; % upper colorado
forc_out = [prefix num2str(centroid(2), '%3.5f') '_' num2str(centroid(1), '%3.5f')];

% directory where to store the lumped inputs
lumped_inputs_dir = '/home/jschap/Documents/ESSD/lumped_cal/'; 

%% Forcings

forcnames = dir([fullfile(forcdir, prefix) '*']);
forc1 = dlmread(fullfile(forcdir, forcnames(1).name), ' ', 0, 0);
ncells = length(forcnames);
nvars = size(forc1, 2);
nt = size(forc1, 1);

[lon, lat] = get_coordinates_from_VIC_file(forcdir, prefix, 'fluxes');

inmask_ind = zeros(ncells,1);
for k=1:ncells
    
    % get indices of cells that are inside the basin mask

    i = find(basin.lon == lon(k));
    j = find(basin.lat == lat(k));
    ind = sub2ind(size(basin.mask), j, i);
    skipflag = logical(basin.mask(ind)); % 0 if not in mask, 1 if in mask
    if skipflag == 0
        inmask_ind(k) = 0;
    else
        inmask_ind(k) = 1;
    end
end

inmask_inds = find(inmask_ind);

forc = zeros(nt, nvars);
count = 0;
for k=inmask_inds'
    newforc = dlmread(fullfile(forcdir, forcnames(k).name), ' ', 1, 0);
    forc = forc + newforc;
    count = count + 1;
    disp(count)    
end
forc_avg = forc/sum(inmask_ind);
save(fullfile(lumped_inputs_dir, 'average_forc.mat'),'forc_avg');

% count = 0;
% for k=1:ncells
%     
%     % skip cells that are outside the basin mask
% 
%     i = find(basin.lon == lon(k));
%     j = find(basin.lat == lat(k));
%     ind = sub2ind(size(basin.mask), j, i);
%     skipflag = logical(basin.mask(ind));
%     if skipflag == 0
%         continue
%     else
%         newforc = dlmread(fullfile(forcdir, forcnames(k).name), ' ', 1, 0);
%         forc = forc + newforc;
%         count = count + 1;
%         disp(count)
%     end
% end
% forc_avg = forc/count;
% save(fullfile(lumped_inputs_dir, 'average_forc.mat'),'forc_avg');

% Write out the lumped forcings
fID_FORC = fopen(fullfile(lumped_inputs_dir, forc_out), 'w');
% formatspec = '%3.5f %3.5f %3.5f %3.5f %3.5f %3.5f %3.5f\n';
formatspec = '%3.5f %3.5f %3.5f %3.5f %3.5f %3.5f %3.5f %3.5f\n';
fprintf(fID_FORC, formatspec, forc_avg');
fclose(fID_FORC);

%% Soils

soils_lumped = mean(basin.soils);
soils_lumped(1) = 1; % run cell
soils_lumped(2) = 1; % grid cell
soils_lumped(3) = centroid(2); % lat
soils_lumped(4) = centroid(1); % lon
soils_lumped(53) = 1; % fs_active
write_soils(5, soils_lumped, fullfile(lumped_inputs_dir, 'soils_ucrb_lumped.txt'), '3l');

%% Elevation bands

elevation_bands = load('/home/jschap/Documents/Data/VICGlobal/Classic/snowbands_MERIT.txt');

% Find the elevation bands associated with the grid cells in the basin
[aa, bb] = ismember(elevation_bands(:,1), basin.soils(:,2));

% Keep those bands only
basin.elevation_bands = elevation_bands(aa, :);

% Take their average
elevbands_lumped = mean(basin.elevation_bands);
elevbands_lumped(:,1) = 1; % gridcell ID

% Write out the new (lumped) elevation band file
elevbands_outname = fullfile(lumped_inputs_dir, 'elev_bands_ucrb_lumped.txt');
dlmwrite(elevbands_outname, elevbands_lumped, 'Delimiter', '\t','precision', 8);

%% Vegetation parameters

% This one can take some time because it has to load the vegetation
% parameter file

% vegfile = '/home/jschap/Documents/Data/VICGlobal/Classic/global_vegetation_1_16_IGBP.txt';
% ncells = size(soils, 1);
% savename = '/Volumes/HD3/VICParametersGlobal/Global_1_16/vegetation/Vegetation_Fractions/mats/vegparamtable_IGBP_11252019.mat';
% Note: read_vegparam saves nvegtable, etc. under savename
% [nvegtable, vegparamtable, latlontable, LC] = read_vegparam(vegfile, soils, ncells, savename);

A = load('/home/jschap/Documents/Data/VICGlobal/Data/Vegetation/vegparamtable_IGBP.mat');
nvegtable = A.nvegtable;
vegparamtable = A.vegparamtable;

[bb, ~] = ismember(nvegtable(:,1), basin.soils(:,2)); % this is identical to aa, above

vpt_lumped = struct();
names1 = fieldnames(vegparamtable);
nvegclasses = length(names1);
vegtype_exists = zeros(nvegclasses, 1);
cover_fraction = zeros(nvegclasses, 1);
for i=1:nvegclasses
    vpt_lumped.(names1{i}) = vegparamtable.Water_Bodies(1,:);
    
    params_for_current_lc_class = vegparamtable.(names1{i})(bb,:);
    
    % eliminate empty rows from consideration
    elim_ind = params_for_current_lc_class(:,2) == 0;
    params_for_current_lc_class(:,1) = 1;
    params_for_current_lc_class = params_for_current_lc_class(~elim_ind,:);
    vpt_lumped.(names1{i}) = mean(params_for_current_lc_class,1);
    
    cover_fraction(i) = vpt_lumped.(names1{i})(2);
    rd1(i) = vpt_lumped.(names1{i})(3);
    rf1(i) = vpt_lumped.(names1{i})(4);
    rd2(i) = vpt_lumped.(names1{i})(5);
    rf2(i) = vpt_lumped.(names1{i})(6);
    rd3(i) = vpt_lumped.(names1{i})(7);
    rf3(i) = vpt_lumped.(names1{i})(8);
    
    if ~isnan(vpt_lumped.(names1{i})(1))
        vegtype_exists(i) = 1;
    end
    
end
nveg = sum(vegtype_exists);

%% Write vegetation parameter file

savename = fullfile(lumped_inputs_dir, 'vegpar_ucrb_lumped.txt');
cellID = 1;

vgn = [17,1,2,3,4,5,6,7,9,8,10,11,12,13,14,15,16];

cover_fraction = cover_fraction(logical(vegtype_exists));
rd1 = rd1(logical(vegtype_exists));
rd2 = rd2(logical(vegtype_exists));
rd3 = rd3(logical(vegtype_exists));
rf1 = rf1(logical(vegtype_exists));
rf2 = rf2(logical(vegtype_exists));
rf3 = rf3(logical(vegtype_exists));
vgn = vgn(logical(vegtype_exists));

% normalize cover fraction to one
cv_norm = cover_fraction/sum(cover_fraction);

% normalize root fractions to one
rf_sum = rf1' + rf2' + rf3';
rf1_norm = rf1'./rf_sum;
rf2_norm = rf2'./rf_sum;
rf3_norm = rf3'./rf_sum;

fID = fopen(savename, 'w');
fmt = '%d %d\n';
fprintf(fID, fmt, [cellID, nveg]);

% For each vegetation class, write parameters
for ii=1:nveg
    current_vgn = vgn(ii);
    fmt = '%d %4.3f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f\n';
    vegpars = [current_vgn cv_norm(ii), rd1(ii), rf1_norm(ii), rd2(ii), rf2_norm(ii), rd3(ii), rf3_norm(ii)];
    fprintf(fID, fmt, vegpars);
end

fclose(fID);

%% Glaciers

% Get basin average glacier cover

[glaciers_pandoh, R_pandoh] = geotiffread('/Volumes/HD3/SWOTDA/Calibration/Pandoh/pandoh_glacier_fract.tif');
glacier_fract = sum(glaciers_pandoh(:))/ncells;
save('/Volumes/HD3/SWOTDA/Calibration/Lumped/Pandoh/pandoh_glacier_fract.mat', 'glacier_fract')
