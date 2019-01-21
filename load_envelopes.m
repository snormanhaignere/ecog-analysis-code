function [D, t, stim_names, electrode_research_numbers] = load_envelopes(exp, subjid, ...
    analysis_name, freq_cutoffs, resp_win, varargin)

% 2017-01-26: Modified to make averaging across repetitions within a run optional
% 
% 2019-01-21: Modified to make baseline period used to calculate PSC values
% a parameter

global root_directory;

%% Parameters

% whether or to overwrite the existing results
I.overwrite = false;

% whether or not to plot results
I.plot = false;

% runs to use
I.runs = [];

% whether or not to average stimulus repetitions within runs
I.average_reps = true;

% debug mode
I.keyboard = false;

% convert to units of percent signal change
I.convert_to_psc = true;

% time window before stimulus onset used to calculate PSC values
I.basewin = 0.5;

% optionally overwrite defaults
I = parse_optInputs_keyvalue(varargin, I);

if I.keyboard
    keyboard;
end

%% Misc setup

% directory for this project
project_directory = [root_directory '/' exp];

% string identifying the frequency range
freq_string = [num2str(freq_cutoffs(1)) '-' num2str(freq_cutoffs(2)) 'Hz'];

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
    
    input_MAT_file = [project_directory '/analysis/envelopes/' ...
        subjid '/r' num2str(I.runs(i)) '/envelopes_' freq_string '_' ...
        analysis_name '_mapped2stim_win' ...
        num2str(resp_win(1)) 'to' num2str(resp_win(2)) '.mat'];
    load(input_MAT_file, 'envelopes_mapped_to_stim', 'env_sr', ...
        'outliers_mapped_to_stim', 'stim_names', 'electrode_research_numbers');
    
    % average across any repetitions with a run
    if i == 1
        [n_smps, n_stimuli, ~, n_electrodes] = size(envelopes_mapped_to_stim);
        if I.average_reps
            D = nan([n_smps, n_stimuli, n_runs, n_electrodes]);
        else
            D = nan([n_smps, n_stimuli, 0, n_electrodes]);
        end
    end
    assert(length(electrode_research_numbers)==n_electrodes);
    
    keyboard;
    envelopes_mapped_to_stim(outliers_mapped_to_stim==1) = NaN;
    
    if I.average_reps
        D(:,:,i,:) = nanmean(envelopes_mapped_to_stim,3);
    else
        D = cat(3, D, envelopes_mapped_to_stim);
    end
    
end

t = resp_win(1):1/env_sr:resp_win(2);

%% Convert to PSC

if I.convert_to_psc
    
    % convert to % signal change by subtracting and dividing by baseline
    % average across all samples and stimuli
    assert(length(t) == n_smps);
    baseline_average = nanmean(nanmean(D(t < 0 & t >= -I.basewin,:,:,:),1),2);
    baseline_average_rep = repmat(baseline_average, [n_smps, n_stimuli, 1, 1]);
    D = 100 * (D - baseline_average_rep) ./ baseline_average_rep;
    
end