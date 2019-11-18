% Set up parallel
%
% Subsets the VIC soil parameter file into multiple, smaller soil parameter files 
% and generates corresponding copies of the VIC global parameter file for running 
% VIC 4 or VIC 5 Classic in parallel as a job array on Hoffman2. 
% This works because computations for one grid cell are independent of other grid cells.
%
% Written by Dongyue Li 4/2019
% Modified 5/1/2019 JRS
% Updated 11/18/2019 JRS to be a function, not a script

function set_up_parallel(control_params)

%% Inputs

% global_param_file = '/Volumes/HD4/SWOTDA/Data/IRB/VIC/FullDomain/globalparam_WB_SB_1980-2018.txt';
% soil_parameter_file = '/Volumes/HD4/SWOTDA/Data/IRB/VIC/FullDomain/soils_snapped.SB';

savedir = fullfile(fileparts(control_params.vic_out_dir), 'parallel');
savedir = '/Volumes/HD4/SWOTDA/Data/IRB/VIC/FullDomain/vic_parallel_WB_1980-2018'; % directory to save the newly created global param files
mkdir(savedir)
savedir_soil = '/u/home/j/jschaper/vic/IRB/vic_parallel'; % directory to save the subsetted soil param files
savedir_out = '/u/flashscratch/j/jschaper/Output/'; % directory where VIC results should be output

soil_basename = './soils_irb'; % basename for subdivided soil parameter files
global_basename = './global_param_irb';

ncells = 75231;
nl = 1000; % number of grid cells to keep in each soil parameter file
n = ceil(ncells/nl); % how many divisions to divide the soil parameter file into

% Create n directories so the outputs from each division are saved in different directories (for debugging purposes)
outdir_names = cell(n,1);
for jj=1:n
  outdir_names{jj} = fullfile(savedir_out, ['Output_' num2str(jj)]);
%   mkdir(outdir_names{jj})
end

%% Make copies of the global parameter file

for jj=1:n

    % read the global parameter file
    fid = fopen(global_param_file, 'r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);

    % change the soil parameter line
    A{153} = ['SOIL ' fullfile(savedir_soil, [soil_basename, '_', num2str(jj), '.txt'])];

    % change the output directory
    A{182} = ['RESULT_DIR ' outdir_names{jj}];

    % Write cell A into txt
    sn = fullfile(savedir, [global_basename, '_', num2str(jj), '.txt']);
    
    fid = fopen(sn, 'w');
    for i=1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);

end

%% Subdivide the soil parameter file

aa=dlmread(soil_parameter_file);
ncells = size(aa,1);
jj=76;
for jj=1:n

    bb=[];
     
    if jj==n
        bb = aa((jj-1)*nl+1:ncells,:);
    else
        bb=aa((jj-1)*nl+1:(jj-1)*nl+nl,:);
    end

    fid = fopen(fullfile(savedir, [soil_basename, '_', num2str(jj) '.txt']), 'w');

    for i = 1:size(bb,1)

        % May need to modify the format specs here, depending on the variables included in the soil parameter file

        fprintf(fid,'%i %i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i %i %i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i %i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %i %i %i %i %i %i\n', bb(i,:));

    end

    fclose(fid);

end
