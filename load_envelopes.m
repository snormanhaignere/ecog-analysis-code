function [D, t, stim_names, electrode_research_numbers] = ...
    load_envelopes(exp, subjid, preproc_idstring, varargin)

% 2017-01-26: Modified to make averaging across repetitions within a run optional
% 
% 2019-01-21: Modified to make baseline period used to calculate PSC values
% a parameter, made it possible to resample/rewindow as well as interpolate
% NaN value

global root_directory;

%% Parameters

% whether or to overwrite the existing results
I.overwrite = false;

% whether or not to plot results
I.plot = false;

% runs to use
I.runs = [];

% whether or not to average stimulus repetitions within runs
I.avreps = false;

% debug mode
I.keyboard = false;

% convert to units of percent signal change
I.psc = true;

% time window before stimulus onset used to calculate PSC values
I.basewin = 0.5;

% can optionally re-window
I.win = [];

% can optionally resample to new rate
I.sr = NaN;

% whether or not to interpolate NaN values
I.interpNaNs = true;

% optionally overwrite defaults
I = parse_optInputs_keyvalue(varargin, I);

if I.keyboard
    keyboard;
end

%% Misc setup

% directory for this project
project_directory = [root_directory '/' exp];

%% Load data

% load the envelope data
% average reps within runs, concatenate across runs
% D, outiers:
% -> samples/time x stimuli x runs x electrodes
if isempty(I.runs)
    I.runs = read_runs(exp, subjid);
end
n_runs = length(I.runs);
for i = 1:n_runs
    
    fprintf('Loading run %d\n', i);
    drawnow;
    
    % load MAT file with preprocessed envelopes
    input_MAT_file = [project_directory '/analysis/preprocessing/' ...
            subjid '/r' num2str(I.runs(i)) '/' preproc_idstring '/env.mat'];
    load(input_MAT_file, 'envelopes_mapped_to_stim', 'env_sr', ...
        'outliers_mapped_to_stim', 'stim_names', 'electrode_research_numbers', 'resp_win');
    
    % average across any repetitions with a run
    if i == 1
        [n_smps, n_stimuli, ~, n_electrodes] = size(envelopes_mapped_to_stim);
        if I.avreps
            D = nan([n_smps, n_stimuli, n_runs, n_electrodes]);
        else
            D = nan([n_smps, n_stimuli, 0, n_electrodes]);
        end
    end
    assert(length(electrode_research_numbers)==n_electrodes);
    
    envelopes_mapped_to_stim(outliers_mapped_to_stim>0.1) = NaN;
    
    if I.avreps
        D(:,:,i,:) = nanmean(envelopes_mapped_to_stim,3);
    else
        D = cat(3, D, envelopes_mapped_to_stim);
    end
    
end

t = resp_win(1):1/env_sr:resp_win(2);

%% Interpolate NaN values

% interpolate NaN values
if I.interpNaNs
    D = interpNaN_ndarray(1:n_smps, D, 'pchip');
end

%% Convert to PSC

if I.psc
    
    % convert to % signal change by subtracting and dividing by baseline
    % average across all samples and stimuli
    assert(length(t) == n_smps);
    baseline_average = nanmean(nanmean(D(t < 0 & t >= -I.basewin,:,:,:),1),2);
    baseline_average_rep = repmat(baseline_average, [n_smps, n_stimuli, 1, 1]);
    D = 100 * (D - baseline_average_rep) ./ baseline_average_rep;
    
end

%% Resample/rewindow

% resample to desired window and sample rate
if ~isempty(I.win) && any(abs(I.win-resp_win)>1e-3) ...
    || ~isnan(I.sr) && abs(I.sr - env_sr) > 1e-3
    if isempty(I.win); I.win = resp_win; end
    if isnan(I.sr); I.sr = env_sr; end
    D = resample_and_window(D, resp_win, env_sr, I.win, I.sr);
    t = I.win(1):1/I.sr:I.win(2);
end