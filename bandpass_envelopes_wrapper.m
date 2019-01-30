function [env_MAT_file, param_idstring] = bandpass_envelopes_wrapper(exp, subjid, r, preproc_idstring, varargin)

% Calculates envelopes of bandpassed ECoG signals.
%
% signal_MAT_file is a cell array with one MAT file per run. The MAT file must
% contain a matrix 'signal' with the input signal, a variable 'sr' with the
% signal's sampling rate.
%
% 2016-08-14: Created, Sam NH
%
% 2016-09-23: Minor changes, Sam NH
%
% 2016-10-18: Removed hashing, added analysis name, Sam NH
%
% 2016-10-18: Preprocessing parameter structure now an argument, Sam NH
% 
% 2019-01-21: Changed parameter handling

%% Setup

I.cutoffs = [70, 140];
I.envsr = 100;
I.order = 6;
I.overwrite = false;
[I, ~, C_value] = parse_optInputs_keyvalue(varargin, I);

% string to identify the parameters of envelope analysis
env_idstring = ['env_' struct2string(I,'include_fields', {'cutoffs', 'envsr'}) ...
    '_' struct2string(C_value, 'omit_field', {'overwrite', 'cutoffs', 'envsr'})];
if env_idstring(end)=='_'; env_idstring = env_idstring(1:end-1); end

global root_directory

% directory for this project
project_directory = [root_directory '/' exp];

% combine preprocessing and envelope idstring
param_idstring = [preproc_idstring '/' env_idstring];

% directory to save results to
input_directory = [project_directory '/analysis/preprocessing/' ...
    subjid '/r' num2str(r) '/' preproc_idstring];
output_directory = [project_directory '/analysis/preprocessing/' ...
    subjid '/r' num2str(r) '/' param_idstring];
if ~exist(output_directory, 'dir'); mkdir(output_directory); end

% directory to save figures to
figure_directory = strrep(output_directory, 'analysis', 'figures');
if ~exist(figure_directory, 'dir'); mkdir(figure_directory); end

%% loop through runs

% number of filters
env_MAT_file = [output_directory '/env.mat'];

% check if mat file already exists
if ~exist(env_MAT_file, 'file') || I.overwrite
    
    % load signal matrix
    preproc_MAT_file = [input_directory '/cleaned.mat'];
    load(preproc_MAT_file, 'signal', 'sr', 'good_channels', 'electrode_research_numbers');
    
    % measure envelopes
    fprintf('Extracting envelopes for band %s...\n', ...
        [num2str(I.cutoffs(1)) '-' num2str(I.cutoffs(2)) 'Hz']);
    envelopes = bandpass_envelopes(signal, sr, I.envsr, I.cutoffs, I.order,...
        figure_directory, 'good_channels', good_channels, ...
        'electrode_numbers', electrode_research_numbers); %#ok<*NASGU>
    
    % save to file
    env_sr = I.envsr;
    cutoffs = I.cutoffs;
    save(env_MAT_file, 'envelopes', 'env_sr', 'cutoffs', ...
        'good_channels', 'electrode_research_numbers', 'r', '-v7.3');
    
end
