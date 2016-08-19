function [gamma_stimulus_response, gamma_stimulus_response_t] = map_gamma_envelopes_to_stimulus_onsets(gamma_env, gamma_env_t, stim_onsets_sec, response_window)

% function [gamma_stimulus_response, gamma_stimulus_response_t] = map_gamma_envelopes_to_stimulus_onsets(gamma_env, gamma_env_t, stim_onsets_sec, response_window)
% 
% Subdivides a vector of gamma envelope responses to many stimuli into chunked/stimulus locked
% responses to all of the stimuli present.
% 
% -- Inputs --
% 
% gamma_env: vector of gamma envelope magnitudes
% 
% env_sr: sampling rate of the gamma envelopes
% 
% stim_onsets_sec: vector of stimulus onsets in seconds
% 
% response_window: range of times relative to stimulus onset to extract
% 
% -- Outputs -- 
% 
% gamma_stimulus_response: [samples/time x stimulus-onsets x electrode] matrix of gamma response
% timecourses in response to each stimulus
% 
% gamma_stimulus_response_t: corresponding time stamps, equal to:
% (response_window(1):1/env_sr:response_window(2))'
% 
% 2016-1-26: Created by Sam NH

% envelope sampling rate
env_sr = 1/(gamma_env_t(2)-gamma_env_t(1));

% number of stimulus onsets
n_stimulus_onsets = length(stim_onsets_sec);

% timestamps for each gamma power measurement within the response window
gamma_stimulus_response_t = (response_window(1):1/env_sr:response_window(2))';

% number of samples within the window
n_smps = length(gamma_stimulus_response_t);

% number of electrodes
n_electrodes = size(gamma_env,2);

% sample x stimulus-onset x electrode
target_times = repmat(gamma_stimulus_response_t, 1, n_stimulus_onsets) + repmat(stim_onsets_sec', n_smps, 1);
gamma_stimulus_response_unwrapped = interp1(gamma_env_t, gamma_env, target_times(:), 'pchip');
gamma_stimulus_response = reshape(gamma_stimulus_response_unwrapped, [n_smps, n_stimulus_onsets, n_electrodes]);