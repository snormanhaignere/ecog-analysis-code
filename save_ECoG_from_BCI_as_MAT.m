function MAT_file = save_ECoG_from_BCI_as_MAT(exp, subjid, r, varargin)

% Saves ECoG data from BCI .dat file as a .mat file
%
% 2016-09-23: Created, Sam NH
% 
% 2018-03-03: Added debug mode, make it possible to skip trigger assignment

global root_directory;

I.overwrite = false;
I.electrode_order = [];
I.keyboard = false;
I.skip_audio_trigger = false;
I = parse_optInputs_keyvalue(varargin, I);

% debug mode
if I.keyboard
    keyboard;
end

% directory for this project
project_directory = [root_directory '/' exp];

% directory with the data for this experiment
data_directory = [project_directory '/data/ECoG/' subjid];

% directory to save results to
analysis_directory = [project_directory '/analysis/preprocessing/' subjid '/r' num2str(r)];
if ~exist(analysis_directory, 'dir')
    mkdir(analysis_directory);
end

% directory to save figures to
figure_directory = strrep(analysis_directory,'analysis','figures');
if ~exist(figure_directory, 'dir')
    mkdir(figure_directory);
end

% check if mat file already exists
MAT_file = [analysis_directory '/raw.mat'];
if ~exist(MAT_file, 'file') || I.overwrite
    
    % load the raw data and parameters
    % convert to double precision
    bci_run_file = [data_directory '/r' num2str(r) '.dat'];
    fprintf('Loading signal...\n'); drawnow;
    [signal, states, parameters, total_samples] = load_bcidat(bci_run_file);
    
    % convert to double
    signal = double(signal);
    
    % optionally chance ordering of electrodes
    if ~isempty(I.electrode_order)
        signal = signal(:,I.electrode_order);
    end
    
    % save sampling rate as separate variable
    sr = parameters.SamplingRate.NumericValue;
    
    % separate out and save audio trigger
    if ~I.skip_audio_trigger
        load([data_directory '/electrode_types.mat'], 'audio_trigger');
        if isnumeric(audio_trigger)
            audio_trigger_signal = signal(:,audio_trigger);
        elseif ischar(audio_trigger)
            audio_trigger_signal = states.(audio_trigger);
        else
            error('Conditional fell through');
        end
    else
        audio_trigger_signal = [];
    end
    
    % select ECoG electrodes
    load([data_directory '/electrode_types.mat'], 'ecog');
    electrode_research_numbers = ecog;
    signal = signal(:,ecog);    

    % save as MAT file
    save(MAT_file, 'signal', 'sr', 'states', ...
        'parameters', 'total_samples', ...
        'electrode_research_numbers', 'audio_trigger_signal');
    
end

