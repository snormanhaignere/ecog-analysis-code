function [MAT_file_with_envelopes_mapped_to_stim, param_idstring] = ...
    map_envelopes_to_stimuli(exp, subjid, r, env_idstring, resp_win, ...
    fixed_duration, varargin)

% Maps a vector of envelope values to the onset of stimuli. Essentially a
% wrapper for map_signal_to_stimulus_onsets.m.
%
% 2016-08-19: Created, Sam NH
%
% 2016-09-23: Minor Changes, Sam NH
% 
% 2018-11-19: Made it possible to have variable duration stimuli
% 
% 2019-01-21: Made fixed/variable duration a parameter that must be
% specified

I.overwrite = false;
I.remove1backs = false;
I.durtol = 1e-2;
I.stimfile = 'stim_names';
[I, ~, C_value] = parse_optInputs_keyvalue(varargin, I);

global root_directory

% directory for this project
project_directory = [root_directory '/' exp];

% string to identify the parameters of this analysis
map2stim_idstring = [...
    'win' num2str(resp_win(1)) 'to' num2str(resp_win(2)) '_fixdur' num2str(fixed_duration) ...
    '_' struct2string(C_value, 'omit_field', {'overwrite', 'durtol'})];
if map2stim_idstring(end) == '_'; map2stim_idstring = map2stim_idstring(1:end-1); end
param_idstring = [env_idstring '/map2stim_' map2stim_idstring];

% directory to save results to
input_directory = [project_directory '/analysis/preprocessing/' subjid '/r' num2str(r) '/' env_idstring];
output_directory = [project_directory '/analysis/preprocessing/' subjid '/r' num2str(r) '/' param_idstring];
if ~exist(output_directory, 'dir'); mkdir(output_directory); end

% directory to save figures to
figure_directory = strrep(output_directory, 'analysis', 'figures');
if ~exist(figure_directory, 'dir'); mkdir(figure_directory); end

% file with mapped envelopes
MAT_file_with_envelopes_mapped_to_stim = [output_directory '/env.mat'];

% check if mat file already exists
if ~exist(MAT_file_with_envelopes_mapped_to_stim, 'file') || I.overwrite
    
    % names of all of the stimuli in the experiment
    load([root_directory '/' exp '/analysis/' I.stimfile '.mat'],'stim_names');
           
    % load the envelopes
    env_MAT_file = [input_directory '/env.mat'];
    load(env_MAT_file, 'envelopes', 'outliers', 'env_sr', ...
        'good_channels', 'electrode_research_numbers');
    
    % information about the onset of stimuli in each run
    % NULL periods are removed
    para_file = [root_directory '/' exp '/data/para/' subjid '/r' num2str(r) '.par'];
    t = stim_timing_from_para(para_file, 'remove-NULL');
    
    % can optionally remove back to back presentations of the same stimulus
    if I.remove1backs
        rep1back = false(1,length(t.names));
        for i = 2:length(t.names)
            rep1back(i) = strcmp(t.names{i-1}, t.names{i});
        end
        t.names(rep1back) = [];
        t.ons_in_sec(rep1back) = [];
        t.dur_in_sec(rep1back) = [];
        t.ids(rep1back) = [];
        clear rep1back;
    end
    
    % map to onsets 
    if fixed_duration
        
        assert(max(abs(t.dur_in_sec-t.dur_in_sec(1))) < I.durtol)
        
        % map the enevelopes to the stimuli
        envelopes_mapped_to_stim = map_signal_to_stimulus_onsets(...
            envelopes, env_sr, stim_names, t, resp_win);
        
        % map the outliers to the stimuli
        outliers_mapped_to_stim = map_signal_to_stimulus_onsets(...
            double(outliers), env_sr, stim_names, t, resp_win);
        
    else
        
        % assert(max(abs(t.dur_in_sec-t.dur_in_sec(1))) > I.durtol)
        
        % map the enevelopes to the stimuli
        envelopes_mapped_to_stim = map_signal_to_stimulus_onsets_variable_duration(...
            envelopes, env_sr, stim_names, t, resp_win);
        
        % map the outliers to the stimuli
        outliers_mapped_to_stim = map_signal_to_stimulus_onsets_variable_duration(...
            double(outliers), env_sr, stim_names, t, resp_win);

    end
    
    % save results
    save(MAT_file_with_envelopes_mapped_to_stim, ...
        'envelopes_mapped_to_stim', 'outliers_mapped_to_stim',...
        'env_sr', 'resp_win', 'stim_names', 'good_channels', ...
        'electrode_research_numbers', '-v7.3');
    
end
