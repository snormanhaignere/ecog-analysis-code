function mapped_signal = map_signal_to_stimulus_onsets(...
    signal, sr, stim_names, stim_onsets_sec, response_window)

% function [mapped_signal, t] = map_signal_to_stimulus_onsets(...
%     signal, sr, stim_onsets_sec, response_window)
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
% stim_onsets_sec: vector of stimulus onsets in seconds
% 
% response_window: range of times relative to stimulus onset to extract
% 
% -- Outputs -- 
% 
% mapped_signal: [samples/time x stimulus-onsets x electrode] matrix of signal
% values
% 
% 2016-1-26: Created by Sam NH

assert(ndims(signal)<=2);
assert(isvector(stim_onsets_sec));
assert(isscalar(sr));
assert(length(response_window)==2);

% ensure column vector
stim_onsets_sec = stim_onsets_sec(:);

% timing vector
% time x stimulus
t = (response_window(1):1/sr:response_window(2))';
n_t = length(t);
n_stimulus_onsets = length(stim_onsets_sec);
target_times = repmat(t, 1, n_stimulus_onsets) ...
    + repmat(stim_onsets_sec', n_t, 1);

% interpolate times
% (sample * stimulus onset) x electrode
mapped_signal_unwrapped = ...
    interp1((0:length(signal)-1)'/sr, signal, target_times(:), 'pchip');

% reshape so as to separate out different stimuli
n_electrodes = size(signal,2);
mapped_signal = reshape(mapped_signal_unwrapped, ...
    [n_t, n_stimulus_onsets, n_electrodes]);