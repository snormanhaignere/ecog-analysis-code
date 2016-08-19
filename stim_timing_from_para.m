function t = stim_timing_from_para(para_file, varargin)

% Reads the onset (t.ons_in_sec), duration (t.dur_in_sec), stimulus name
% (t.names), and stimulus id (t.ids) of each event in a paradigm file.
% 
% 2016-08-19 - Created, Sam NH

% write to para file
fid = fopen(para_file, 'r');

% read file contents
file_contents = textscan(fid,'%f%d%f%f%s');
fclose(fid);
[t.ons_in_sec, t.ids, t.dur_in_sec, ~, t.names] = file_contents{:}; %#ok<*NASGU>

% optionally remove NULL periods
if optInputs(varargin, 'remove-NULL')
    xi = ~ismember(t.names, 'NULL');
    t.ons_in_sec = t.ons_in_sec(xi);
    t.dur_in_sec = t.dur_in_sec(xi);
    t.names = t.names(xi);
    t.ids = t.ids(xi);
end


