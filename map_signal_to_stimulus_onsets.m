function signal_mapped_by_stim = map_signal_to_stimulus_onsets(...
    signal, sr, stim_names, S, response_window)

% function [signal_mapped_by_stim, t] = map_signal_to_stimulus_onsets(...
%     signal, sr, stim_onsets, response_window)
%
% Subdivides a signal vector so that the responses are locked to the onset
% of a stimulus / task.
%
% -- Inputs --
%
% signal: [signal x electrode] matrix of envelopes
%
% sr: signal sampling rate
%
% stim_names: cell array with the name for each stimulus presented in the
% experiment
%
% S: structure containing the time stamp of each stimulus onset
% (S.ons_in_sec), the duration (S.dur_in_sec), and the
% corresponding stimulus name (S.name); the duration is not used
%
% response_window: range of times relative to stimulus onset to extract
%
% -- Outputs --
%
% signal_mapped_by_stim: [samples/time x stim x repetition x electrode] matrix of signal
% values
%x
% 2016-1-26: Created by Sam NH
%
% 2016-08-19 - Elaborated so that stimuli with the same name are grouped
% together, Sam NH

assert(ndims(signal)<=2);
assert(isvector(S.ons_in_sec));
assert(isscalar(sr));
assert(length(response_window)==2);

% ensure column vector
S.ons_in_sec = S.ons_in_sec(:);
n_onsets = length(S.ons_in_sec);

% timing vector
% time x stimulus
t = (response_window(1):1/sr:response_window(2))';
n_t = length(t);
target_times = repmat(t, 1, n_onsets) ...
    + repmat(S.ons_in_sec', n_t, 1);

% interpolate times
% (sample * stimulus onset) x electrode
signal_mapped_by_onset = ...
    interp1((0:length(signal)-1)'/sr, signal, target_times(:), 'pchip');

% reshape so as to separate out different stimuli
% -> sample x stim onset x electrode
n_electrodes = size(signal,2); 
signal_mapped_by_onset = reshape(signal_mapped_by_onset, ...
    [n_t, n_onsets, n_electrodes]);

% group responses with the same stimulus name
% determine number of reps per stimulus
n_stimuli = length(stim_names);
stim_index_for_each_onset = nan(1, n_onsets);
n_reps_per_stim = zeros(1, n_stimuli);
for i = 1:n_onsets
    if any(strcmp(S.names{i}, stim_names))
        stim_index_for_each_onset(i) = find(strcmp(S.names{i}, stim_names));
        n_reps_per_stim(stim_index_for_each_onset(i)) = ...
            n_reps_per_stim(stim_index_for_each_onset(i)) + 1;
    end
end
% assert(all(~isnan(stim_index_for_each_onset)));

% create the final matrix
% -> sample x stim x rep x electrode
n_electrodes = size(signal,2);
signal_mapped_by_stim = nan(n_t, n_stimuli, max(n_reps_per_stim), n_electrodes);
rep = zeros(1, n_stimuli);
for i = 1:n_onsets
    if any(strcmp(S.names{i}, stim_names))
        rep(stim_index_for_each_onset(i)) = rep(stim_index_for_each_onset(i))+1;
        signal_mapped_by_stim(:, stim_index_for_each_onset(i), rep(stim_index_for_each_onset(i)), :) ...
            = signal_mapped_by_onset(:,i,:);
    end
end