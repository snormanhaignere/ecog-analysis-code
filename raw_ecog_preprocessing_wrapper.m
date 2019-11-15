function [MAT_file_with_preproc_signal, param_idstring] = ...
    raw_ecog_preprocessing_wrapper(exp, subjid, r, varargin)

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

I.steps = {'60Hz', 'car', 'notch'};
I.bw60Hz = 0.6;
I.frac = 0.5; % see channel_selection_from_60Hz_noise.m
I.min_nchannels = 10; % see channel_selection_from_60Hz_noise.m
I.thresh60Hz = 5;
I.hpcutoff = 0.5;
I.hporder = 4;
I.notchbw = 1;
I.notchfreqs = [60, 120, 180, 240];
I.exclude_from_car = [];
I.array_inds = [];
I.overwrite = false;
[I, ~, C_value] = parse_optInputs_keyvalue(varargin, I);

% string to identify the parameters of this analysis
param_idstring = struct2string(C_value, 'omit_field', {'overwrite','array_inds'});
if ~isempty(I.array_inds); param_idstring = [param_idstring '_arrayspecified']; end
if isempty(param_idstring); param_idstring = 'default'; end

% directory for this project
project_directory = [root_directory '/' exp];

% directory to save results to
input_directory = [project_directory '/analysis/preprocessing/' subjid '/r' num2str(r) '/'];
output_directory = [input_directory '/' param_idstring];
if ~exist(output_directory, 'dir'); mkdir(output_directory); end

% directory to save figures to
figure_directory = strrep(output_directory,'analysis','figures');
if ~exist(figure_directory, 'dir'); mkdir(figure_directory); end

% check if mat file already exists
MAT_file_with_preproc_signal = [output_directory '/cleaned.mat'];
if ~exist(MAT_file_with_preproc_signal, 'file') || I.overwrite
    
    % load the raw data and sampling rate
    fprintf('Loading signal...\n'); drawnow;
    load([input_directory '/raw.mat'], ...
        'signal', 'sr', 'electrode_research_numbers');
    
    % preprocess the data    
    [signal, good_channels, noise60Hz_rms] = ...
        raw_ecog_preprocessing(signal, sr, figure_directory, ...
        'steps', I.steps, ...
        'bw60Hz', I.bw60Hz, 'frac', I.frac, 'thresh60Hz', I.thresh60Hz, ...
        'hpcutoff', I.hpcutoff, 'hporder', I.hporder, ...
        'notchbw', I.notchbw, 'notchfreqs', I.notchfreqs, ...
        'exclude_from_car', I.exclude_from_car, ...
        'electrode_numbers', electrode_research_numbers, ...
        'array_inds', I.array_inds, ...
        'min_nchannels', I.min_nchannels);
    
    % save to mat file
    save(MAT_file_with_preproc_signal, ...
        'signal', 'electrode_research_numbers', 'good_channels', 'noise60Hz_rms', 'sr', 'r', '-v7.3');
end



