function [ons_in_sec, dur_in_sec, stim_names, stim_ids] ...
    = stim_onsets_from_bci(bci_states, bci_parameters)

% function [ons_in_sec, dur_in_sec, stim_names, stim_ids] =
% stim_onsets_from_bci(bci_states, bci_parameters)
% 
% Calculates stimulus onsets and corresponding stimulus names from the states
% and parameters structures of a bci file (returned by the function
% load_bcidat.m)
% 
% -- Inputs --
% 
% bci_states, bci_parameters: structures returned by load_bcidat function, see
% example below
% 
% -- Outputs --
% 
% stim_onsets_in_sec: onsets of each stimulus in seconds
% 
% stim_names: corresponding names for each stimulus
% 
% stim_ids: a unique number assigned to each stimulus in BCI
% 
% -- Example -- bci_run_file = '/Users/svnh2/Dropbox
% (MIT)/mcdexp-svnh/ecog-analysis-code/example_dataset.dat'; [~, bci.states,
% bci.parameters, ~] = load_bcidat(bci_run_file); [stim_onsets_sec, stim_names,
% stim_ids] = stim_onsets_from_bci(bci.states, bci.parameters);
% 
% 2016-08-19 - Modified to measure durations of each period, NULL periods are no
% longer removed but flagged as such

bci_states.StimulusCode = double(bci_states.StimulusCode);

% detect onsets in samples
onsets_logical = [bci_states.StimulusCode(1)~=0; diff(bci_states.StimulusCode)~=0];
onsets_smps = find(onsets_logical);
durations_smps = diff([onsets_smps; length(bci_states.StimulusCode)]);

% corresponding ids at the onset of each stimulus
stim_ids = bci_states.StimulusCode(onsets_smps);

% % remove null onsets
% xi = stim_ids ~= 0;
% onsets_smps = onsets_smps(xi);
% durations_smps = durations_smps(xi);
% stim_ids = stim_ids(xi);

% corresponding names for each id
stim_names = cell(1, length(stim_ids));
stim_names(stim_ids~=0) = bci_parameters.Stimuli.Value(1,stim_ids(stim_ids~=0));
stim_names(stim_ids==0) = {'NULL'};

% convert samples to seconds
raw_ecog_sr = bci_parameters.SamplingRate.NumericValue;
ons_in_sec = (onsets_smps-1)/raw_ecog_sr;
dur_in_sec = durations_smps/raw_ecog_sr;





