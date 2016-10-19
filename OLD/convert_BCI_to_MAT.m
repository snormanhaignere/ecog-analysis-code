function MAT_file = convert_BCI_to_MAT(exp, subjid, r, varargin)

% Converts a BCI file to a MAT file
% 
% 2016-08-12 - Created, Sam NH

% general-purpose ecog analysis code
global root_directory;

% directory for this project
project_directory = [root_directory '/' exp];

% directory with the data for this experiment
data_directory = [project_directory '/data/ECoG/' subjid];

% directory to save results to
analysis_directory = [project_directory '/analysis/preprocessing/' subjid];
if ~exist(analysis_directory, 'dir');
    mkdir(analysis_directory);
end

% directory to save figures to
figure_directory = strrep(analysis_directory,'analysis','figures');
if ~exist(figure_directory, 'dir');
    mkdir(figure_directory);
end

% check if mat file already exists
MAT_file = [data_directory '/r' num2str(r) '.mat'];
if exist(MAT_file, 'file') && ~optInputs(varargin, 'overwrite')
    return;
end

% load the raw data and parameters
% convert to double precision
bci_run_file = [data_directory '/r' num2str(r) '.dat'];
fprintf('Loading signal...\n'); drawnow;
[signal, states, parameters, total_samples, file_samples ] ...
    = load_bcidat(bci_run_file); %#ok<ASGLU>
signal = double(signal);
sr = parameters.SamplingRate.NumericValue;

% select out the ECoG electrodes
load([data_directory '/electrode_types.dat'], 'ecog');
electrode_research_numbers = ecog;
signal = signal(:,ecog);

% save as MAT file
save(MAT_file, 'signal', 'sr', 'states', ...
    'parameters', 'total_samples', 'file_samples', 'electrode_research_numbers');


