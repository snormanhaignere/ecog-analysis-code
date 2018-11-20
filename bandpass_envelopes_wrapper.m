function MAT_file_with_envelopes = bandpass_envelopes_wrapper(...
    exp, subjid, r, analysis_name, P, varargin)

% Calculates envelopes of bandpassed ECoG signals.
%
% signal_MAT_file is a cell array with one MAT file per run. The MAT file must
% contain a matrix 'signal' with the input signal, a variable 'sr' with the
% signal's sampling rate.
%
% 2016-08-14: Created, Sam NH
%
% 2016-09-23: Minor changes, Sam NH
% 
% 2016-10-18: Removed hashing, added analysis name, Sam NH
% 
% 2016-10-18: Preprocessing parameter structure now an argument, Sam NH

%% Setup

I.overwrite = false;
I = parse_optInputs_keyvalue(varargin, I);

global root_directory

% directory for this project
project_directory = [root_directory '/' exp];

% directory to save results to
analysis_directory = [project_directory '/analysis/envelopes/' subjid '/r' num2str(r)];
if ~exist(analysis_directory, 'dir')
    mkdir(analysis_directory);
end

% directory to save figures to
figure_directory = strrep(analysis_directory, 'analysis', 'figures');
if ~exist(figure_directory, 'dir')
    mkdir(figure_directory);
end

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
        '/envelopes_' bp_freq_range_string '_' analysis_name '.mat'];
    
    % check if mat file already exists
    if ~exist(MAT_file_with_envelopes{i}, 'file') || I.overwrite
        
        % load signal matrix
        if ~signal_matrix_loaded
            signal_MAT_file = [project_directory '/analysis/preprocessing' ...
                '/' subjid '/r' num2str(r) '/cleaned_signal_' analysis_name '.mat'];
            load(signal_MAT_file, 'signal', 'sr', 'good_channels');
        end
        
        % measure envelopes
        fprintf('Extracting envelopes for band %s...\n', bp_freq_range_string);
        envelopes = bandpass_envelopes(signal, sr, P.bandpass_env_sr, ...
            P.bandpass_cutoffs_in_Hz(:,i), P.bandpass_filter_orders(i),...
            figure_directory); %#ok<*NASGU>
        
        % save to file
        env_sr = P.bandpass_env_sr;
        band_in_Hz = P.bandpass_cutoffs_in_Hz(:,i);
        save(MAT_file_with_envelopes{i}, 'envelopes', 'env_sr', 'band_in_Hz', ...
            'good_channels', 'P', 'r', '-v7.3');
        
    end
    
end
