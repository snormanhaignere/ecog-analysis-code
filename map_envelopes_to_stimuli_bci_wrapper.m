function MAT_file_with_envelopes_mapped_to_stimuli = ...
    map_envelopes_to_stimuli_bci_wrapper(exp, subjid, envelopes_MAT_file, varargin)

%% Setup

clc;

global root_directory

% directory for this project
project_directory = [root_directory '/' exp];

% directory with the data for this experiment
data_directory = [project_directory '/data/ECoG/' subjid '/'];

% directory to save results to
analysis_directory = [project_directory '/analysis' ...
    '/envelopes-mapped-to-stim/' subjid];
if ~exist(analysis_directory, 'dir');
    mkdir(analysis_directory);
end

% names for all of the stimuli tested
load([project_directory '/code/naturalsound_stim_names.mat'],'stim_names');
n_stimuli = 165;

for i = 1:n_runs
    
    load(envelopes_MAT_file, 'r');
    
    matfile = [analysis_directory '/r' num2str(r) '.mat'];
    if exist(matfile, 'file') && ~optInputs(varargin, 'overwrite')
        return;
    end
    
    %%
    
    fprintf('Loading run parameters...\n')
    [~,bci.states,bci.parameters,~] = ...
        load_bcidat([data_directory '/r' num2str(r) '.dat']);
    
    % onsets of stimuli in seconds
    [stim_onsets_sec, R.stim_names, stim_ids] = ...
        stim_onsets_from_bci(bci.states, bci.parameters);
    
    % remove repeated stimulus presentations
    repeated_blocks = find([false; diff(stim_ids)==0]);
    R.stim_names(repeated_blocks) = [];
    stim_onsets_sec(repeated_blocks) = [];
    stim_ids(repeated_blocks) = []; %#ok<NASGU>
    
    response_window = [-1 3];
    [R.gamma_stimulus_response, gamma_stimulus_response_t] = ...
        map_gamma_envelopes_to_stimulus_onsets(gamma_env, gamma_env_t, stim_onsets_sec, response_window);
    
    %%
    
    fprintf('Run %d...\n',runs(i)); drawnow;
    
    fprintf('Loading preprocessed data...\n');
    load([project_directory '/analysis/raw-ecog-preprocessing/'...
        schalk_subjid '/r' num2str(runs(i)) '.mat'],...
        'preproc_signal', 'ecog_sr', 'good_channels');
    

    
    %% Map to envelope for each stimulus
    
    
    % initialize based on first run
    if i == 1
        [n_smps, n_stim_onsets, n_electrodes] = size(R.gamma_stimulus_response);
        gamma_stimulus_response = nan([n_smps, n_stimuli, n_runs, n_electrodes]);
    end
    
    % map each stimulus to a fixed location in the matrix
    for j = 1:n_stimuli
        xi = ismember(R.stim_names, stim_names{j}); %#ok<USENS>
        if sum(xi) ~= 1
            error('Should be exactly one stimulus per run');
        end
        gamma_stimulus_response(:,j,i,:) = R.gamma_stimulus_response(:,xi,:);
    end
    
end

% save results
save(matfile, 'gamma_stimulus_response', ...
    'gamma_stimulus_response_t', 'stim_names', ...
    'env_sr', 'ecog_sr', 'good_channels');