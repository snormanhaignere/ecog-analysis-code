function runs = read_runs(exp, subjid, varargin)

% 2016-08-14: Created, Sam NH
% 
% 2019-11-11: Modifed to use para directory

% general-purpose ecog analysis code
global root_directory;

% directory for this project
project_directory = [root_directory '/' exp];

% directory with the data for this experiment
para_directory = [project_directory '/data/para/' subjid '/'];

% parse run numbers from files in data directory
files_in_data_directory = mydir(para_directory, '.par');
runs = [];
for i = 1:length(files_in_data_directory)
    runstr = regexp(files_in_data_directory{i}, 'r(\d)+.par', 'match');
    if ~isempty(runstr)
        runs = [runs, str2double(runstr{1}(2:end-4))]; %#ok<AGROW>
    end
end

% % directory with the data for this experiment
% data_directory = [project_directory '/data/ECoG/' subjid '/'];
% 
% % parse run numbers from files in data directory
% files_in_data_directory = mydir(data_directory, '.dat', '_aux');
% if isempty(files_in_data_directory)
%     files_in_data_directory = mydir(data_directory, '.mat', '_aux');
% end
% runs = [];
% for i = 1:length(files_in_data_directory)
%     runstr = regexp(files_in_data_directory{i}, 'r(\d)+', 'match');
%     if ~isempty(runstr)
%         runs = [runs, str2double(runstr{1}(2:end))]; %#ok<AGROW>
%     end
% end