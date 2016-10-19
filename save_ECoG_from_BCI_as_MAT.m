function MAT_file = save_ECoG_from_BCI_as_MAT(exp, subjid, r, varargin)

% Saves ECoG data from BCI .dat file as a .mat file
%
% 2016-09-23: Created, Sam NH

global root_directory;

I.overwrite = false;
I.electrode_order = [];
I = parse_optInputs_keyvalue(varargin, I);

% directory for this project
project_directory = [root_directory '/' exp];

% directory with the data for this experiment
data_directory = [project_directory '/data/ECoG/' subjid];

% directory to save results to
analysis_directory = [project_directory '/analysis/preprocessing/' subjid '/r' num2str(r)];
if ~exist(analysis_directory, 'dir');
    mkdir(analysis_directory);
end

% directory to save figures to
figure_directory = strrep(analysis_directory,'analysis','figures');
if ~exist(figure_directory, 'dir');
    mkdir(figure_directory);
end

% check if mat file already exists
MAT_file = [analysis_directory '/raw.mat'];
if ~exist(MAT_file, 'file') || I.overwrite
    
    % load the raw data and parameters
    % convert to double precision
    bci_run_file = [data_directory '/r' num2str(r) '.dat'];
    fprintf('Loading signal...\n'); drawnow;
    [ signal, states, parameters, total_samples, file_samples ] ...
        = load_bcidat(bci_run_file); %#ok<ASGLU>
    
    % convert to double
    signal = double(signal);
    
    % optionally chance ordering of electrodes
    if ~isempty(I.electrode_order)
        signal = signal(:,I.electrode_order);
    end
    
    % save sampling rate as separate variable
    sr = parameters.SamplingRate.NumericValue; %#ok<NASGU>
    
    % select ECoG electrodes
    load([data_directory '/electrode_types.mat'], 'ecog');
    electrode_research_numbers = ecog; %#ok<NASGU>
    signal = signal(:,ecog); %#ok<NASGU>
    
    % save as MAT file
    save(MAT_file, 'signal', 'sr', 'states', ...
        'parameters', 'total_samples', ...
        'file_samples', 'electrode_research_numbers');
    
end

