function runs = read_runs(exp, subjid, varargin)

% 2016-08-14 - Created, Sam NH 

% general-purpose ecog analysis code
global root_directory;

% directory for this project
project_directory = [root_directory '/' exp];

% directory with the data for this experiment
data_directory = [project_directory '/data/ECoG/' subjid '/'];

% parse run numbers from files in data directory
files_in_data_directory = mydir(data_directory);
runs = [];
for i = 1:length(files_in_data_directory)
    runstr = regexp(files_in_data_directory{i}, 'r(\d)+\.dat', 'match');
    if ~isempty(runstr)
        runs = [runs, str2double(regexp(runstr{1}, '(\d)+', 'match'))]; %#ok<AGROW>
    end
end