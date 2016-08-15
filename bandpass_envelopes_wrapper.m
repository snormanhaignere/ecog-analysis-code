function MAT_file_with_envelopes = bandpass_envelopes_wrapper(exp, subjid, r, signal_MAT_file, varargin)

% Calculates envelopes of bandpassed ECoG signals.
%
% signal_MAT_file is a cell array with one MAT file per run. The MAT file must
% contain a matrix 'signal' with the input signal, a variable 'sr' with the
% signal's sampling rate.
%
% 2016-08-14 - Created, Sam NH

%% Setup

global root_directory

% directory for this project
project_directory = [root_directory '/' exp];

% directory to save results to
analysis_directory = [project_directory '/analysis/preprocessing/' subjid];
if ~exist(analysis_directory, 'dir');
    mkdir(analysis_directory);
end

% directory to save figures to
figure_directory = strrep(analysis_directory, 'analysis', 'figures');
if ~exist(figure_directory, 'dir');
    mkdir(analysis_directory);
end

% parameters of the filters
P = preprocessing_parameters;

% create a hash string specific to the inputs and parameters to this function
all_args = [exp, subjid, r, signal_MAT_file, varargin];
assert(length(all_args) == nargin);
relevant_parameters = {P.bandpass_env_sr, P.bandpass_cutoffs_in_Hz, ...
    P.bandpass_filter_orders};
hash_string = DataHash([all_args, relevant_parameters]);

%% loop through runs

% number of filters
n_bands = size(P.bandpass_cutoffs_in_Hz,2);
signal_matrix_loaded = false;
MAT_file_with_envelopes = cell(1, n_bands);
for i = 1:n_bands
    
    % MAT file to save results to
    bp_freq_range_string = ...
        [num2str(P.bandpass_cutoffs_in_Hz(1,i)) ...
        '-' num2str(P.bandpass_cutoffs_in_Hz(2,i)) 'Hz'];
    MAT_file_with_envelopes{i} = [analysis_directory ...
        '/r' num2str(r) '_bpfilt_' bp_freq_range_string '_' hash_string '.mat'];
    
    % check if mat file already exists
    if exist(MAT_file_with_envelopes{i}, 'file') && ~optInputs(varargin, 'overwrite')
        continue;
    end
    
    % load signal matrix
    if ~signal_matrix_loaded
        load(signal_MAT_file, 'signal', 'sr', 'good_channels');
    end
    
    % measure envelopes
    fprintf('Extracting envelopes for band %s...\n', bp_freq_range_string);
    envelopes = bandpass_envelopes(signal, sr, P.bandpass_env_sr, ...
        P.bandpass_cutoffs_in_Hz(:,i), P.bandpass_filter_orders(i),...
        figure_directory, ['r' num2str(r) '_bpfilt_' bp_freq_range_string], ...
        good_channels, varargin); %#ok<*NASGU>
    
    % save to file
    env_sr = P.bandpass_env_sr;
    band_in_Hz = P.bandpass_cutoffs_in_Hz(:,i);
    save(MAT_file_with_envelopes{i}, 'envelopes', 'env_sr', 'band_in_Hz', ...
        'good_channels', 'P', 'r');
    
end
