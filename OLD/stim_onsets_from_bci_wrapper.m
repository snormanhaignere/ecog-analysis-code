function stim_onsets = stim_onsets_from_bci_wrapper(exp, subjid, r, varargin)

% function stim_onsets = stim_onsets_from_bci_wrapper(exp, subjid, r, varargin)
% 
% Wrapper for stim_onsets_from_bci that determines the onsets for each run of
% this experiment.
% 
% 2016-08-19 - Created Sam NH 

global root_directory

% directory for this project
project_directory = [root_directory '/' exp];

% directory with the data for this experiment
data_directory = [project_directory '/data/ECoG/' subjid '/'];

% load parameters
[~,bci.states,bci.parameters,~] = ...
    load_bcidat([data_directory '/r' num2str(r) '.dat']);

% determine onsets
[stim_onsets.time_sec, stim_onsets.names, stim_onsets.ids] = ...
    stim_onsets_from_bci(bci.states, bci.parameters);
