% Read bounds for soil parameters
%
% INPUTS
% Column indices of soil parameter file for calibration parameters
% Lookup table with upper and lower bounds
%
% OUTPUTS
% Upper and lower bounds for calibration parameters

function [bl, bu] = read_bounds(params, parbounds)

fID = fopen(parbounds, 'r');
dat = textscan(fID, '%s%s%s%s%s', 'delimiter', '\t', 'headerlines', 1);
fclose(fID);

bl = str2double(dat{3}(params))';
bu = str2double(dat{4}(params))';

return