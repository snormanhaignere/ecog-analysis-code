function MAT_file_with_preproc_signal = ...
    raw_ecog_preprocessing_wrapper(exp, subjid, r, varargin)

% Preprocessing scripts applied to the raw ecog data stored as a MAT file.
% Acts as a wrapper / file-handler for raw_ecog_preprocessing.m. Returns mat
% files with the preprocessed ECoG signals for each run. As
%
% 2016-08-12: Created, Sam NH
%
% 2016-09-23: Changed how optional inputs are handled, Sam NH

global root_directory;

I.overwrite = false;
I = parse_optInputs_keyvalue(varargin, I);

% parameters of the analysis
P = preprocessing_parameters;

% directory for this project
project_directory = [root_directory '/' exp];

% directory to save results to
analysis_directory = [project_directory '/analysis/preprocessing/' subjid '/r' num2str(r)];
if ~exist(analysis_directory, 'dir');
    mkdir(analysis_directory);
end

% directory to save figures to
figure_directory = strrep(analysis_directory,'analysis','figures');
if ~exist(figure_directory, 'dir');
    mkdir(figure_directory);
end

% create a hash string specific to the inputs and parameters to this function
hash_string = DataHash({exp, subjid, r, P.bw_60Hz_peak_filt, P.hp_filt_order, ...
    P.hp_filt_cutoff_in_Hz, P.notch_n_harmonics, P.notch_bw});

% check if mat file already exists
MAT_file_with_preproc_signal = [analysis_directory ...
    '/cleaned_signal_' hash_string '.mat'];
if ~exist(MAT_file_with_preproc_signal, 'file') || I.overwrite
    
    % load the raw data and sampling rate
    fprintf('Loading signal...\n'); drawnow;
    load([analysis_directory '/raw.mat'], ...
        'signal', 'sr', 'electrode_research_numbers');
    
    % preprocess the data
    [signal, good_channels, noise60Hz_rms] = ...
        raw_ecog_preprocessing(signal, sr, P, figure_directory,...
        'electrode_numbers', electrode_research_numbers); %#ok<NODEF,ASGLU>
    
    % save to mat file
    save(MAT_file_with_preproc_signal, ...
        'signal', 'good_channels', 'noise60Hz_rms', 'sr', 'r');
    
end



