% Reads a parameter file, skipping blank lines and comments
%
% Sometime, adapt this to a function for general use reading
% ASCII parameter files

function B = read_parfile(fname)

% Read the parameter file (from read_global_paramfile.m)
fid = fopen(fname, 'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

% Get the non-comment, non-space lines
nlines = length(A);
ind = zeros(nlines, 1);
for k=1:nlines
    if length(A{k})>=1
        if ~strcmp(A{k}(1), '#')
            ind(k) = 1;
        end
    end
end
nlines_w_values = sum(ind);
B = cell(nlines_w_values, 1);
count1 = 1;
for k=1:nlines
    if ind(k)==1
        B{count1} = A{k};
        count1 = count1+1;
    end
end

B = B(1:(nlines_w_values-1));

return

%% For future consideration

% regexp(B{1}, '\w+\s/\w{20}*', 'match')
% 
% dirname = regexp(B{1}, '(\w+)\s/(\w+/?\w+)', 'tokens')
% 
% 
% filecontent = fileread(parameter_file);
% 
% filenames = '/*\.\w\w\w';
% filenames = '^\w+\s\/\w+\.\w+$';
% filenames = '/\w+?/\w';
% files = regexp(filecontent, filenames, 'match')
% files = regexp(filecontent, '(?<=/.)[^.]*$', 'match')
% 
% directories = regexp(filecontent, '^\s(/\w+\w*)', 'tokens', 'lineanchors')
% 
% comments = regexp(filecontent, '^#*', 'match', 'lineanchors')
% 
% paramvalues = regexp(filecontent, '^\t*(\w+)\s*=\s*(\w+)\s*$', 'tokens', 'lineanchors');
% paramvalues = vertcat(paramvalues{:});
% [path1, name1, ext1] = fileparts(ss{2});

% Specialized application (assumes order of pars in parfile)