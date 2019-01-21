function envelope_outliers_wrapper(exp, subjid, r, envelope_MAT_file, varargin)

% Detects outliers in ECoG envelopes. A wrapper for envelope_outliers.m. See
% this function for details.
%
% 2016-08-15 - Created, Sam NH
% 
% 2018-11-20 - Saves the multiscale outliers, includes optional parameters

I.overwrite = false;
I.scales = [0 0.25 1 4];
I.thresholds = [5 2.5 2 1.5];
I.percentile = 0.9;
I.plot_mosaic = true;
I.plot_individ_electrodes = false;
I = parse_optInputs_keyvalue(varargin, I);

% general-purpose ecog analysis code
global root_directory;

% directory for this project
project_directory = [root_directory '/' exp];

% check if file already exists
outliers_exist = ~isempty(whos('-file', envelope_MAT_file, 'outliers'));
if ~outliers_exist || I.overwrite
    
    % load envelopes and envelope sampling rate
    load(envelope_MAT_file, 'envelopes', 'env_sr', 'band_in_Hz');
    
    % MAT file to save results to
    bp_freq_range_string = ...
        [num2str(band_in_Hz(1)) '-' num2str(band_in_Hz(2)) 'Hz'];
    
    % directory to save plots to
    figure_directory = [project_directory '/figures/envelopes/' subjid '/' ...
        'r' num2str(r) '/bpfilt-' bp_freq_range_string '-outliers'];
    if ~exist(figure_directory, 'dir')
        mkdir(figure_directory);
    end
    
    % detect outliers
    fprintf('Detecting outliers...\n'); drawnow;
    [outliers, multiscale_outliers] = envelope_outliers(envelopes, ...
        env_sr, figure_directory, 'scales', I.scales, 'thresholds', I.thresholds, ...
        'percentile', I.percentile, 'plot_mosaic', I.plot_mosaic, ...
        'plot_individ_electrodes', I.plot_individ_electrodes);
    
    % save results
    save(envelope_MAT_file, 'outliers', 'multiscale_outliers', '-append', '-v7.3');
    
end