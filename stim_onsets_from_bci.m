function [stim_onsets_sec, stim_names, stim_ids] = stim_onsets_from_bci(bci_states, bci_parameters)

% function [stim_onsets_sec, stim_ids] = stim_onsets_from_bci(bci_states, bci_parameters, raw_ecog_sr)
% 
% Calculates stimulus onsets and corresponding stimulus names from the states and parameters
% structures of a bci file (returned by the function load_bcidat.m)
% 
% -- Inputs --
% 
% bci_states, bci_parameters: structures returned by load_bcidat function, see example below
% 
% -- Outputs --
% 
% stim_onsets_sec: onsets of each stimulus in seconds
% 
% stim_names: corresponding names for each stimulus
% 
% stim_ids: a unique number assigned to each stimulus in BCI
% 
% -- Example --
% bci_run_file = '/Users/svnh2/Dropbox (MIT)/mcdexp-svnh/ecog-analysis-code/example_dataset.dat';
% [~, bci.states, bci.parameters, ~] = load_bcidat(bci_run_file);
% [stim_onsets_sec, stim_names, stim_ids] = stim_onsets_from_bci(bci.states, bci.parameters);

% detect onsets in samples
stimulus_onsets_logical = [bci_states.StimulusCode(1)~=0; diff(double(bci_states.StimulusCode))~=0];
stimulus_onsets_smps = find(stimulus_onsets_logical);

% corresponding ids at the onset of each stimulus
stim_ids = double(bci_states.StimulusCode(stimulus_onsets_smps));

% remove null onsets
xi = stim_ids ~= 0;
stimulus_onsets_smps = stimulus_onsets_smps(xi);
stim_ids = stim_ids(xi);

% corresponding names for each id
stim_names = bci_parameters.Stimuli.Value(1,stim_ids);

% convert samples to seconds
raw_ecog_sr = bci_parameters.SamplingRate.NumericValue;
stim_onsets_sec = (stimulus_onsets_smps-1)/raw_ecog_sr;




