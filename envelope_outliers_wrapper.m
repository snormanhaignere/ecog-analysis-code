function envelope_outliers_wrapper(exp, subjid, r, envelope_MAT_file, varargin)

% Detects outliers in ECoG envelopes. A wrapper for envelope_outliers.m. See
% this function for details.
% 
% 2016-08-15 - Created, Sam NH

% general-purpose ecog analysis code
global root_directory;

% directory for this project
project_directory = [root_directory '/' exp];

% directory to save plots to
figure_directory = [project_directory '/figures/preprocessing/' subjid];

% number of standard deviations above the mean to exclude
P = preprocessing_parameters;

% check if file already exists
outliers_exist = ~isempty(whos('-file', envelope_MAT_file, 'outliers'));
if outliers_exist && ~optInputs(varargin, 'overwrite')
    return
end

% load envelopes and envelope sampling rate
load(envelope_MAT_file, 'envelopes');

% detect outliers
fprintf('Detecting outliers...\n'); drawnow;
outliers = envelope_outliers(envelopes, ...
    P.envelope_outlier_threshold, figure_directory, ['r' num2str(r)]); %#ok<NASGU>

% save results
save(envelope_MAT_file, 'outliers', '-append');

