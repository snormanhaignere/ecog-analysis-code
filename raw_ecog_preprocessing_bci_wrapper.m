function MAT_file_with_preproc_signal = ...
    raw_ecog_preprocessing_bci_wrapper(exp, subjid, r, varargin)

% Preprocessing scripts applied to the raw ecog data stored as a bcidat file.
% Acts as a wrapper / file-handler for raw_ecog_preprocessing.m. Returns mat
% files with the preprocessed ECoG signals for each run. As
%
% 2016-08-12 - Created, Sam NH

% preprocessing parameters
P = preprocessing_parameters;

% general-purpose ecog analysis code
global root_directory;

% directory for this project
project_directory = [root_directory '/' exp];

% directory with the data for this experiment
data_directory = [project_directory '/data/ECoG/' subjid '/'];

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
MAT_file_with_preproc_signal = [analysis_directory ...
    '/r' num2str(r) '_cleaned_signal.mat'];
if exist(MAT_file_with_preproc_signal, 'file') && ~optInputs(varargin, 'overwrite')
    return;
end

% load the raw data and parameters
% convert to double precision
bci_run_file = [data_directory '/r' num2str(r) '.dat'];
fprintf('Loading signal...\n'); drawnow;
[signal, ~, params, ~] = load_bcidat(bci_run_file);
signal = double(signal);
sr = params.SamplingRate.NumericValue;

% preprocess the data
[signal, good_channels, noise60Hz_rms] = ...
    raw_ecog_preprocessing(signal, sr, P, ...
    figure_directory, ['r' num2str(r)], varargin{:}); %#ok<ASGLU>

% save to mat file
save(MAT_file_with_preproc_signal, ...
    'signal', 'good_channels', 'noise60Hz_rms', 'sr', 'r');




