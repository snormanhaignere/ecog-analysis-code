function para_file = save_stim_timings_in_bci_as_para(exp, subjid, r, varargin)

% function para_file = save_stim_timings_in_bci_as_para(exp, subjid, r, varargin)
%
% Extracts the stimulus timing info contained in the BCI files and saves this
% information as a paradigm file.

% 2016-08-19 - Created Sam NH
%
% 2018-03-02: Modified to make it optional as to whether the stimulus names are
% read from the parameter structure in the BCI file. The alternative is to read
% the names from the analysis/stim_names.mat, in which case the order of the
% stimuli in the MAT file should match the indices.
%
% 2018-03-21: Modified so that one can add both constant and stimulus-specific
% delays to the onset of each stimulus
% 
% 2018-11-20: Previously the stim_id in the para file could potentially
% change. To fix this the index is now tied to the index of each stimulus
% name in the stim_names.mat file (or the corresponding delay file).

global root_directory

I.fn_to_stim_name = @(x)x;
I.overwrite = false;
I.stim_names_from_bci = true;
I.constant_delay = 0;
I.stim_spec_delays = false;
I.keyboard = false;
I = parse_optInputs_keyvalue(varargin, I);

if I.keyboard
    keyboard;
end

% directory for this project
project_directory = [root_directory '/' exp];

% directory with the data for this experiment
data_directory = [project_directory '/data/ECoG/' subjid '/'];

% directory to save para files to
para_directory = [project_directory '/data/para/' subjid '/'];
if ~exist(para_directory, 'dir')
    mkdir(para_directory);
end

% para file to write
para_file = [para_directory '/r' num2str(r) '.par'];
if ~exist(para_file, 'file') || I.overwrite
    
    % load parameters
    [~,bci.states,bci.parameters,~] = ...
        load_bcidat([data_directory '/r' num2str(r) '.dat']);
    
    % determine onsets and durations from bci file
    [ons_in_sec, dur_in_sec, stim_name_for_each_onset] = ...
        stim_onsets_from_bci(bci.states, bci.parameters, ...
        'stim_names_from_bci', I.stim_names_from_bci, 'exp', exp);
    
    % load stimulus specific onsets
    if I.stim_spec_delays
        load([project_directory '/analysis/stim_delays.mat'], 'delays', 'stim_names');
        X = load([root_directory '/' exp '/analysis/stim_names.mat'],'stim_names');
        assert(isequal(X.stim_names, stim_names));
    else
        load([root_directory '/' exp '/analysis/stim_names.mat'],'stim_names');
    end
    
    % write to para file
    fid = fopen(para_file, 'w');
    n_onsets = length(ons_in_sec);
    for i = 1:n_onsets
        
        if ~strcmp(I.fn_to_stim_name(stim_name_for_each_onset{i}), 'NULL')
            stim_index = find(ismember(stim_names, I.fn_to_stim_name(stim_name_for_each_onset{i})));
            if isempty(stim_index)
                continue;
            else
                assert(length(stim_index)==1);
            end
        else
            stim_index = 0;
        end
        
        % stimulus specific delay
        if ~strcmp(I.fn_to_stim_name(stim_name_for_each_onset{i}), 'NULL') && I.stim_spec_delays
            stim_spec_delay = delays(stim_index);
        else
            stim_spec_delay = 0;
        end
                
        % write to file
        fprintf(fid,'%10.6f%5d%10.6f%5d%50s\n', ...
            ons_in_sec(i) + I.constant_delay + stim_spec_delay, ...
            stim_index, dur_in_sec(i), NaN, ...
            I.fn_to_stim_name(stim_name_for_each_onset{i}));
        
    end
    fclose(fid);
    
end



