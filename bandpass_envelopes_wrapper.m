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

%% loop through runs

% check if mat file already exists
MAT_file_with_envelopes = [analysis_directory '/r' num2str(r) '.mat'];
if exist(MAT_file_with_envelopes, 'file') && ~optInputs(varargin, 'overwrite')
    return;
end

% measure envelopes
fprintf('Extracting envelopes...')
load(signal_MAT_file, 'signal', 'sr', 'good_channels');
envelopes = bandpass_envelopes(signal, sr, P.bandpass_env_sr, ...
    P.bandpass_cutoffs_in_Hz, P.bandpass_filter_orders,...
    figure_directory, ['r' num2str(r)], good_channels, varargin); %#ok<*NASGU>

% save to file
env_sr = P.bandpass_env_sr;
bands_in_Hz = P.bandpass_cutoffs_in_Hz;
save(MAT_file_with_envelopes, 'envelopes', 'env_sr', 'bands_in_Hz', ...
    'good_channels', 'P', 'r');