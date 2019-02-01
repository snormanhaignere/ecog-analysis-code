function [D, t, stim_names, electrode_research_numbers] = load_envelopes_variable_duration(...
    exp, subjid, preproc_idstring, varargin)

% 2018-11-19: Created, Sam NH

global root_directory;

% whether or to overwrite the existing results
I.overwrite = false;

% whether or not to plot results
I.plot = false;

% runs to use
I.runs = [];

% whether or not to average stimulus repetitions within runs
I.avreps = true;

% debug mode
I.keyboard = false;

% whether or not to interpolate outliers
I.interpNaNs = true;

% convert to units of percent signal change
I.psc = true;

% time-period to treat as baseline
I.basewin = [-0.5, 0];

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
    
    % load data
    % load MAT file with preprocessed envelopes
    input_MAT_file = [project_directory '/analysis/preprocessing/' ...
        subjid '/r' num2str(I.runs(i)) '/' preproc_idstring '/env.mat'];
    load(input_MAT_file, 'envelopes_mapped_to_stim', 'env_sr', ...
        'outliers_mapped_to_stim', 'stim_names', 'electrode_research_numbers', 'resp_win');
    n_stimuli = length(stim_names);
    assert(n_stimuli == length(envelopes_mapped_to_stim));
    assert(n_stimuli == length(outliers_mapped_to_stim));
    
    % initialize
    if i == 1
        D = cell(1, n_stimuli);
        for j = 1:n_stimuli
            [n_smps, ~, n_electrodes] = size(envelopes_mapped_to_stim{j});
            if I.avreps
                D{j} = nan([n_smps, n_runs, n_electrodes]);
            else
                D{j} = nan([n_smps, 0, n_electrodes]);
            end
            clear n_smps;
        end
    end
    
    % exclude outliers
    for j = 1:n_stimuli
        envelopes_mapped_to_stim{j}(outliers_mapped_to_stim{j}>0.1) = NaN;
    end
    
    % average across reps
    for j = 1:n_stimuli
        if I.avreps
            D{j}(:,i,:) = nanmean(envelopes_mapped_to_stim{j},2);
        else
            D{j} = cat(2, D{j}, envelopes_mapped_to_stim{j});
        end
    end
end

% Corresponding time vector
t = cell(1, n_stimuli);
for j = 1:n_stimuli
    t{j} = (0:size(D{j},1)-1)/env_sr+resp_win(1);
end

%% Interpolate NaN values

% interpolate NaN values
if I.interpNaNs
    for j = 1:n_stimuli
        D{j} = interpNaN_ndarray(1:size(D{j},1), D{j}, 'pchip');
    end
end

%% Convert to % signal change

if I.psc
    
    assert(resp_win(1)<0);
    
    % timepoint x electrode (across all stimuli and repetitions)
    baseline_tps = [];
    for j = 1:n_stimuli
        xi = t{j}>I.basewin(1) & t{j}<I.basewin(2);
        baseline_tps = [baseline_tps; ...
            reshape(D{j}(xi,:,:), sum(xi)*size(D{j},2), n_electrodes)]; %#ok<AGROW>
    end
    
    % average all timepoints
    baseline_average = reshape(nanmean(baseline_tps,1), [1, 1, n_electrodes]);
    
    % convert to psc
    for j = 1:n_stimuli
        D{j} = 100 * bsxfun(@times, bsxfun(@minus, D{j}, baseline_average), 1./baseline_average);
    end
end
