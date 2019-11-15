function signal_mapped_by_stim = map_signal_to_stimulus_onsets_variable_duration(...
    signal, sr, stim_names, S, response_window)

% function [signal_mapped_by_stim, t] = map_signal_to_stimulus_onsets_variable_duration(...
%     signal, sr, stim_onsets, response_window)
%
% Subdivides a signal vector so that the responses are locked to the onset
% of a stimulus / task. Similar to map_signal_to_stimulus_onsets, but
% does not assume that each stimulus is the same duration. As a
% consequence, the response window is relative to stimulus onset and
% offset. And there is a different matrix for each stimulus, since the
% matrices can be of different sizes.
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
% response_window: range of times relative to stimulus onset and offset
%
% -- Outputs --
%
% signal_mapped_by_stim: {stim} -> [samples/time x repetition x electrode]
%
% 20189-11-19 - Created, Sam NH

assert(ndims(signal)<=2);
assert(isvector(S.ons_in_sec));
assert(isscalar(sr));
assert(length(response_window)==2);

% ensure column vector
S.ons_in_sec = S.ons_in_sec(:);
n_onsets = length(S.ons_in_sec);

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
if ~all(~isnan(stim_index_for_each_onset));
    keyboard;
end

% assign
n_electrodes = size(signal,2);
rep = zeros(1, n_stimuli);
signal_mapped_by_stim = cell(1, n_stimuli);
for i = 1:n_onsets
    if any(strcmp(S.names{i}, stim_names))

        % repetition index
        rep(stim_index_for_each_onset(i)) = rep(stim_index_for_each_onset(i))+1;

        % interpolate window
        t = (response_window(1):1/sr:S.dur_in_sec(i)+response_window(2))' + S.ons_in_sec(i);
        signal_mapped_by_onset = interp1((0:length(signal)-1)'/sr, signal, t(:), 'pchip');
        
        % initialize
        if rep(stim_index_for_each_onset(i)) == 1
            signal_mapped_by_stim{stim_index_for_each_onset(i)} = nan(length(signal_mapped_by_onset), ...
                n_reps_per_stim(stim_index_for_each_onset(i)), n_electrodes);
        end
        
        % assign
        signal_mapped_by_stim{stim_index_for_each_onset(i)}(:, rep(stim_index_for_each_onset(i)), :) ...
            = signal_mapped_by_onset;
    end
end