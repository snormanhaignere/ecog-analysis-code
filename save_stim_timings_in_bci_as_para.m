function para_file = save_stim_timings_in_bci_as_para(exp, subjid, r, varargin)

% function para_file = save_stim_timings_in_bci_as_para(exp, subjid, r, varargin)
%
% Extracts the stimulus timing info contained in the BCI files and saves this
% information as a paradigm file.
%
% 2016-08-19 - Created Sam NH

global root_directory

I.fn_to_stim_name = @(x)x;
I.overwrite = false;
I = parse_optInputs_keyvalue(varargin, I);

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
    [ons_in_sec, dur_in_sec, stim_name_for_each_onset, stim_ids] = ...
        stim_onsets_from_bci(bci.states, bci.parameters);
    
    % write to para file
    fid = fopen(para_file, 'w');
    n_onsets = length(ons_in_sec);
    for i = 1:n_onsets
        fprintf(fid,'%10.6f%5d%10.6f%5d%50s\n', ...
            ons_in_sec(i), stim_ids(i), dur_in_sec(i), NaN, ...
            I.fn_to_stim_name(stim_name_for_each_onset{i}));
    end
    fclose(fid);
    
end



