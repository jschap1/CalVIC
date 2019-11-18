% Strsplit, but returns a particular part of the result
%
% Specific for use with run_sceua

function result = strsplit_spec(B)

ss = strsplit(B);
result = ss{2};

return