function [gamma_stimulus_response, gamma_stimulus_response_t, stim_names, gamma_env, gamma_env_t] = gamma_stimulus_responses_from_bci(bci_run_file, response_window)

% function [gamma_stimulus_response, gamma_stimulus_response_t, stim_names, gamma_env, gamma_env_t] = gamma_stimulus_responses_from_bci(bci_run_file)
% 
% Wrapper function for other simpler functions that measures stimulus driven gamma power responses
% from a single run of a raw BCI data file. 
% 
% -- Inputs --
% 
% bci_run_file: full path to the bci file
% 
% -- Outputs --
% 
% gamma_stimulus_response: [samples/time x stimulus-onsets x electrode] response matrix matrix containing
% the the raw gamma power response of each electrode to each stimulus onset in the experiment.
% 
% gamma_stimulus_response_t: vector of timestamps for the response matrix above
% 
% stim_names: corresponding stimulus names for the response matrix above
% 
% gamma_env: the full gamma power timecourse, not mapped to any stimulus
% 
% gamma_env_t: corresponding timestamps for the envelope vector above
% 
% 2016-1-26: Created by Sam NH

% load the raw data
[bci.raw_ecog_signal, bci.states, bci.parameters, ~] = load_bcidat(bci_run_file);
raw_ecog_sr = bci.parameters.SamplingRate.NumericValue;

%%

xi = 1:500e3;
stim_time = double(bci.states.StimulusTime(xi));
stim_code = double(bci.states.StimulusCode(xi));
source_time = double(bci.states.SourceTime(xi));

subplot(2,1,1);
plot([stim_time-source_time]);

stim_time = (stim_time - mean(stim_time))/std(stim_time);
stim_code = (stim_code - mean(stim_code))/std(stim_code);
source_time = (source_time - mean(source_time))/std(source_time);

subplot(2,1,2);
plot([stim_time, stim_code, source_time])



%% Stimulus information

% onsets of stimuli in seconds
[stim_onsets_sec, stim_names, stim_ids] = stim_onsets_from_bci(bci.states, bci.parameters);

% remove repeated stimulus presentations
repeated_blocks = find([false; diff(stim_ids)==0]);
stim_names(repeated_blocks) = [];
stim_onsets_sec(repeated_blocks) = [];
stim_ids(repeated_blocks) = []; %#ok<NASGU>

%% Gamma envelopes 

% measure gamma envelopes
env_sr = 100;
gamma_frequency_range = [70 140];
[gamma_env, gamma_env_t] = gamma_envelopes(double(bci.raw_ecog_signal), raw_ecog_sr, env_sr, gamma_frequency_range);

%% Gamma envelopes for each stimulus

[gamma_stimulus_response, gamma_stimulus_response_t] = map_gamma_envelopes_to_stimulus_onsets(gamma_env, gamma_env_t, stim_onsets_sec, response_window);
