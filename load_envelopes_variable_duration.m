function [D, stim_names, env_sr] = load_envelopes_variable_duration(exp, subjid, ...
    analysis_name, freq_cutoffs, resp_win, varargin)

% 2018-11-19: Created, Sam NH

global root_directory;

% window to include when calculating the correlation
I.win = resp_win;

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

% time-period to treat as baseline
I.baseline_period = [-0.5, 0];

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
    
    % load data
    input_MAT_file = [project_directory '/analysis/envelopes/' ...
        subjid '/r' num2str(I.runs(i)) '/envelopes_' freq_string '_' ...
        analysis_name '_mapped2stim_win' ...
        num2str(resp_win(1)) 'to' num2str(resp_win(2)) '.mat'];
    load(input_MAT_file, 'envelopes_mapped_to_stim', 'env_sr', ...
        'outliers_mapped_to_stim', 'stim_names');
    n_stimuli = length(stim_names);
    assert(n_stimuli == length(envelopes_mapped_to_stim));
    assert(n_stimuli == length(outliers_mapped_to_stim));
    
    % initialize
    if i == 1
        D = cell(1, n_stimuli);
        for j = 1:n_stimuli
            [n_smps, ~, n_electrodes] = size(envelopes_mapped_to_stim{j});
            if I.average_reps
                D{j} = nan([n_smps, n_runs, n_electrodes]);
            else
                D{j} = nan([n_smps, 0, n_electrodes]);
            end
        end
    end
    
    % exclude outliers
    for j = 1:n_stimuli
        envelopes_mapped_to_stim{j}(outliers_mapped_to_stim{j}==1) = NaN;
    end
    
    % average across reps
    for j = 1:n_stimuli
        if I.average_reps
            D{j}(:,i,:) = nanmean(envelopes_mapped_to_stim{j},2);
        else
            D{j} = cat(2, D{j}, envelopes_mapped_to_stim{j});
        end
    end
end

%% Convert to % signal change

assert(resp_win(1)<0);

% timepoint x electrode (across all stimuli and repetitions)
baseline_tps = [];
for j = 1:n_stimuli
    t = (0:size(D{j},1)-1)/env_sr+resp_win(1);
    xi = t>I.baseline_period(1) & t<I.baseline_period(2);
    baseline_tps = [baseline_tps; ...
        reshape(D{j}(xi,:,:), sum(xi)*size(D{j},2), n_electrodes)]; %#ok<AGROW>
end

% average all timepoints
baseline_average = reshape(nanmean(baseline_tps,1), [1, 1, n_electrodes]);

% convert to psc
for j = 1:n_stimuli
    D{j} = 100 * bsxfun(@times, bsxfun(@minus, D{j}, baseline_average), 1./baseline_average);
end