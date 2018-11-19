function MAT_file_with_envelopes_mapped_to_stim = ...
    map_envelopes_to_stimuli(exp, subjid, r, envelope_MAT_file, resp_win, varargin)

% Maps a vector of envelope values to the onset of stimuli. Essentially a
% wrapper for map_signal_to_stimulus_onsets.m.
%
% 2016-08-19: Created, Sam NH
%
% 2016-09-23: Minor Changes, Sam NH
% 
% 2018-11-19: Made it possible to have variable duration stimuli

I.overwrite = false;
I.remove_1backs = false;
I.keyboard = false;
I = parse_optInputs_keyvalue(varargin, I);

if I.keyboard
    keyboard;
end

global root_directory;

% file with mapped envelopes
MAT_file_with_envelopes_mapped_to_stim = ...
    strrep(envelope_MAT_file, '.mat', ...
    ['_mapped2stim_win' num2str(resp_win(1)) ...
    'to' num2str(resp_win(2)) '.mat']);

% check if mat file already exists
if ~exist(MAT_file_with_envelopes_mapped_to_stim, 'file') ...
        || I.overwrite
    
    % names of all of the stimuli in the experiment
    load([root_directory '/' exp '/analysis/stim_names.mat'],'stim_names');
    
    % load the envelopes
    load(envelope_MAT_file, 'envelopes', 'outliers', 'env_sr');
    
    % information about the onset of stimuli in each run
    % NULL periods are removed
    para_file = [root_directory '/' exp '/data/para/' subjid '/r' num2str(r) '.par'];
    t = stim_timing_from_para(para_file, 'remove-NULL');
    
    % can optionally remove back to back presentations of the same stimulus
    if I.remove_1backs
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
    if max(abs(t.dur_in_sec-t.dur_in_sec(1))) > 1e-2 % variable duration
        
        % map the enevelopes to the stimuli
        envelopes_mapped_to_stim = map_signal_to_stimulus_onsets_variable_duration(...
            envelopes, env_sr, stim_names, t, resp_win);
        
        % map the outliers to the stimuli
        outliers_mapped_to_stim = map_signal_to_stimulus_onsets_variable_duration(...
            double(outliers), env_sr, stim_names, t, resp_win);
        
    else % fixed duration
        
        % map the enevelopes to the stimuli
        envelopes_mapped_to_stim = map_signal_to_stimulus_onsets(...
            envelopes, env_sr, stim_names, t, resp_win);
        
        % map the outliers to the stimuli
        outliers_mapped_to_stim = map_signal_to_stimulus_onsets(...
            double(outliers), env_sr, stim_names, t, resp_win);
        
    end
    
    % save results
    save(MAT_file_with_envelopes_mapped_to_stim, ...
        'envelopes_mapped_to_stim', 'outliers_mapped_to_stim',...
        'env_sr', 'resp_win', 'stim_names', '-v7.3');
    
end
