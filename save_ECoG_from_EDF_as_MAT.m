function MAT_file = save_ECoG_from_EDF_as_MAT(exp, subjid, r, varargin)

% Saves ECoG data from EDF file as a .mat file
%
% 2019-11-11: Last edited/commented, Sam NH
% 
% 2019-11-11: Updated to optionally exclude a certain window, Sam NH

global root_directory;

I.overwrite = false;
I.ecogchan = [];
I.audiochan = [];
I.trigchan = [];
I.electrode_order = [];
I.keyboard = false;
I.startend = [];
I.excludewin = [];
I.plot = false;
I = parse_optInputs_keyvalue(varargin, I);

% debug mode
if I.keyboard
    keyboard;
end

% directory for this project
project_directory = [root_directory '/' exp];

% directory with the data for this experiment
data_directory = [project_directory '/data/ECoG-EDF/' subjid];

% directory to save results to
analysis_directory = [project_directory '/analysis/preprocessing/' subjid '/r' num2str(r)];
if ~exist(analysis_directory, 'dir')
    mkdir(analysis_directory);
end

% check if mat file already exists
MAT_file = [analysis_directory '/raw.mat'];
if ~exist(MAT_file, 'file') || I.overwrite
    
    % load the raw data and parameters
    edf_file = [data_directory '/r' num2str(r) '.edf'];
    fprintf('Loading signal...\n'); drawnow;
    [hdr, record] = edfread(edf_file);
    
    % optionally chance ordering of electrodes
    if ~isempty(I.electrode_order)
        record = record(I.electrode_order,:);
    end
    
    % select samples
    if ~isempty(I.startend)
        record = record(:, I.startend(1):I.startend(2));
        if ~isempty(I.excludewin)
            I.excludewin = I.excludewin - (I.startend(1)-1);
        end
    end
    
    % save sampling rate as separate variable
    if ~isempty(I.ecogchan)
        sr = hdr.frequency(I.ecogchan(1));
        assert(all(eq_tol(sr, hdr.frequency(I.ecogchan))));
    else
        sr = hdr.frequency(1);
        assert(all(eq_tol(sr, hdr.frequency)));
    end
    
    % separate out audio and trigger
    if ~isempty(I.audiochan)
        audio_signal = record(I.audiochan,:)';
        audio_MAT_file = [project_directory '/data/ECoG-audio/' subjid '/r' num2str(r) '.mat'];
        if I.plot
            figure;
            plot(audio_signal);
            title('Audio');
        end
        save(mkpdir(audio_MAT_file), 'audio_signal', 'sr');
        if ~isempty(I.excludewin)
            excludewin = I.excludewin;
            save(audio_MAT_file, '-append', 'excludewin');
        end
        clear audio_MAT_file audio_signal;
    end
    if ~isempty(I.trigchan)
        trigger_signal = record(I.trigchan,:)';
        trigger_MAT_file = [project_directory '/data/ECoG-trigger/' subjid '/r' num2str(r) '.mat'];
        if I.plot
            figure;
            plot(trigger_signal);
            title('Trigger');
        end
        save(mkpdir(trigger_MAT_file), 'trigger_signal', 'sr');
        if ~isempty(I.excludewin)
            excludewin = I.excludewin;
            save(trigger_MAT_file, '-append', 'excludewin');
        end
        clear trigger_MAT_file trigger_signal;
    end
    
    % pick out ecog
    if ~isempty(I.ecogchan)
        signal = record(I.ecogchan, :)';
        electrode_research_numbers = I.ecogchan;
    else
        signal = record';
        electrode_research_numbers = 1:size(record,1);
    end
    
    % save as MAT file
    save(mkpdir(MAT_file), 'signal', 'sr', 'electrode_research_numbers', '-v7.3');
    if ~isempty(I.excludewin)
        excludewin = I.excludewin;
        save(MAT_file, '-append', 'excludewin');
    end
end

