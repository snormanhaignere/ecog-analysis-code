function [outlier_MAT_file, param_idstring] = envelope_outliers_wrapper(exp, subjid, r, env_idstring, varargin)

% Detects outliers in ECoG envelopes. A wrapper for envelope_outliers.m. See
% this function for details.
%
% 2016-08-15: Created, Sam NH
%
% 2018-11-20: Saves the multiscale outliers, includes optional parameters
% 
% 2019-01-21: Changed parameter handling, Sam NH

I.scales = 0; %[0 0.25 1 4];
I.thresholds = 5; %[5 2.5 2 1.5];
I.percentile = 0.9;
I.overwrite = false;
I.plot_mosaic = true;
I.plot_individ_electrodes = false;
[I, ~, C_value] = parse_optInputs_keyvalue(varargin, I);

global root_directory

% directory for this project
project_directory = [root_directory '/' exp];

% string to identify the parameters of this analysis
outlier_idstring = struct2string(C_value, 'omit_field', ...
    {'overwrite', 'plot_mosaic', 'plot_individ_electrodes'});
if isempty(outlier_idstring); outlier_idstring = 'default'; end

param_idstring = [env_idstring '_out_' outlier_idstring];

% directory to save results to
input_directory = [project_directory '/analysis/preprocessing/' subjid '/r' num2str(r) '/' env_idstring];
output_directory = [project_directory '/analysis/preprocessing/' subjid '/r' num2str(r) '/' param_idstring];
if ~exist(output_directory, 'dir'); mkdir(output_directory); end

% directory to save figures to
figure_directory = strrep(output_directory, 'analysis', 'figures');
if ~exist(figure_directory, 'dir'); mkdir(figure_directory); end

% file to save results to
outlier_MAT_file = [output_directory '/env.mat'];

% check if file already exists
if ~exist(outlier_MAT_file, 'file') || I.overwrite
    
    % load envelopes and envelope sampling rate
    env_MAT_file = [input_directory '/env.mat'];
    load(env_MAT_file, 'envelopes', 'env_sr');
    
    % detect outliers
    fprintf('Detecting outliers...\n'); drawnow;
    [outliers, multiscale_outliers] = envelope_outliers(envelopes, ...
        env_sr, figure_directory, 'scales', I.scales, 'thresholds', I.thresholds, ...
        'percentile', I.percentile, 'plot_mosaic', I.plot_mosaic, ...
        'plot_individ_electrodes', I.plot_individ_electrodes);
    
    % copy over envelopes
    copyfile(env_MAT_file, outlier_MAT_file, 'f');
    
    % append outliers
    save(outlier_MAT_file, 'outliers', 'multiscale_outliers', '-append');
    
end