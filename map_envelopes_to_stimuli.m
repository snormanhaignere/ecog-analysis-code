function MAT_file_with_envelopes_mapped_to_stim = ...
    map_envelopes_to_stimuli(exp, subjid, r, envelope_MAT_file, resp_win, varargin)

% Maps a vector of envelope values to the onset of stimuli. Essentially a
% wrapper for map_signal_to_stimulus_onsets.m. 
% 
% 2016-08-19 - Created, Sam NH

global root_directory;

% file with mapped envelopes
MAT_file_with_envelopes_mapped_to_stim = ...
    strrep(envelope_MAT_file, '.mat', ...
    ['_mapped2stim_win' num2str(resp_win(1)) ...
    'to' num2str(resp_win(2)) '.mat']);

% check if mat file already exists
if exist(MAT_file_with_envelopes_mapped_to_stim, 'file') ...
        && ~optInputs(varargin, 'overwrite')
    return;
end

% names of all of the stimuli in the experiment
load([root_directory '/' exp '/analysis/stim_names.mat'],'stim_names');

% load the envelopes
load(envelope_MAT_file, 'envelopes', 'env_sr');

% information about the onset of stimuli in each run
% NULL periods are removed
para_file = [root_directory '/' exp '/data/para/' subjid '/r' num2str(r) '.par'];
t = stim_timing_from_para(para_file, 'remove-NULL');

% can optionally remove back to back presentations of the same stimulus
if optInputs(varargin, 'remove-1backs')
    repeated_blocks = find([false; diff(t.ids)==0]);
    t.names(repeated_blocks) = [];
    t.ons_in_sec(repeated_blocks) = [];
    t.ids(repeated_blocks) = [];
end

% map the enevelopes to the stimuli
envelopes_mapped_to_stim = map_signal_to_stimulus_onsets(...
    envelopes, env_sr, stim_names, t, resp_win); %#ok<NASGU>

% save results
save(MAT_file_with_envelopes_mapped_to_stim, ...
    'envelopes_mapped_to_stim', 'env_sr', 'resp_win', 'stim_names');

