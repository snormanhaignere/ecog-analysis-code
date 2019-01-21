function MAT_file_with_preproc_signal = ...
    raw_ecog_preprocessing_wrapper(exp, subjid, r, analysis_name, P, varargin)

% Preprocessing scripts applied to the raw ecog data stored as a MAT file.
% Acts as a wrapper / file-handler for raw_ecog_preprocessing.m. Returns mat
% files with the preprocessed ECoG signals for each run. As
%
% 2016-08-12: Created, Sam NH
%
% 2016-09-23: Changed how optional inputs are handled, Sam NH
% 
% 2016-10-18: Removed hashing, added analysis name, Sam NH
% 
% 2018-01-23: Preprocessing parameters are now an input to the function

global root_directory;

I.overwrite = false;
I.exclude_from_car = [];
I.steps = {'60Hz', 'highpass' 'car', 'notch'};
I = parse_optInputs_keyvalue(varargin, I);

% directory for this project
project_directory = [root_directory '/' exp];

% directory to save results to
analysis_directory = [project_directory '/analysis/preprocessing/' subjid '/r' num2str(r)];
if ~exist(analysis_directory, 'dir')
    mkdir(analysis_directory);
end

% directory to save figures to
figure_directory = strrep(analysis_directory,'analysis','figures');
if ~exist(figure_directory, 'dir')
    mkdir(figure_directory);
end

% check if mat file already exists
MAT_file_with_preproc_signal = [analysis_directory ...
    '/cleaned_signal_' analysis_name '.mat'];
if ~exist(MAT_file_with_preproc_signal, 'file') || I.overwrite
    
    % load the raw data and sampling rate
    fprintf('Loading signal...\n'); drawnow;
    load([analysis_directory '/raw.mat'], ...
        'signal', 'sr', 'electrode_research_numbers');
    
    % preprocess the data    
    [signal, good_channels, noise60Hz_rms] = ...
        raw_ecog_preprocessing(signal, sr, P, figure_directory,...
        'electrode_numbers', electrode_research_numbers, ...
        'exclude_from_car', I.exclude_from_car, ...
        'steps', I.steps);
    
    % save to mat file
    save(MAT_file_with_preproc_signal, ...
        'signal', 'electrode_research_numbers', 'good_channels', 'noise60Hz_rms', 'sr', 'r', '-v7.3');
    
end



