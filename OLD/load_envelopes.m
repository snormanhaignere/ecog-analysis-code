function [D,t,sr] = load_envelopes(exp, subjid, analysis_name, freq_cutoffs, resp_win)

global root_directory;  

% gamma envelope matrix
% time x stimuli x runs x electrodes
runs = read_runs(exp, subjid);
n_runs = length(runs);
sr = nan(1,n_runs);
for i = 1:n_runs

    % read in envelopes
    input_MAT_file = [root_directory '/' exp '/analysis/envelopes/' ...
        subjid '/r' num2str(runs(i)) '/envelopes_' ...
        num2str(freq_cutoffs(1)) '-' num2str(freq_cutoffs(2)) 'Hz_' ...
        analysis_name '_mapped2stim_win' ...
        num2str(resp_win(1)) 'to' num2str(resp_win(2)) '.mat'];
    load(input_MAT_file, 'envelopes_mapped_to_stim', 'env_sr', ...
        'outliers_mapped_to_stim', 'stim_names');
    sr(i) = env_sr;
    
    % pre-allocate
    if i == 1
        [n_smps, n_stimuli, ~, n_electrodes] = size(envelopes_mapped_to_stim);
        D = nan([n_smps, n_stimuli, n_runs, n_electrodes]);
    end
    
    % remove outliers
    envelopes_mapped_to_stim(outliers_mapped_to_stim == 1) = NaN; %#ok<AGROW>
    
    % average across reps within a run
    D(:,:,i,:) = nanmean(envelopes_mapped_to_stim,3);
    
end

assert(length(unique(sr))==1);
sr = unique(sr);

% time-vector
t = (0:n_smps-1)/env_sr + resp_win(1);