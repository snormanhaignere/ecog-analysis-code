function outliers_in_stim_mapped_envelopes(exp, subjid, r, envelope_MAT_file, varargin)

% Detects outliers in ECoG envelopes, after the envelopes have been mapped to
% the stimuli. Essentially a wrapper for envelope_outliers.m. See this function
% for details.
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

% load envelopes
load(envelope_MAT_file, 'envelopes_mapped_to_stim');

% unwrap the signal from each electrode
dims = size(envelopes_mapped_to_stim); %#ok<NODEF>
envelopes_mapped_to_stim = reshape(...
    envelopes_mapped_to_stim, [prod(dims(1:3)), dims(4)]);
assert(all(~isnan(envelopes_mapped_to_stim(:))));

% detect outliers
fprintf('Detecting outliers...\n'); drawnow;
outliers = envelope_outliers(envelopes_mapped_to_stim, ...
    P.envelope_outlier_threshold, figure_directory, ...
    ['r' num2str(r) '_mapped_envelopes']);

% reshape to the original dimensions of the matrix
outliers = reshape(outliers, dims); %#ok<NASGU>

% save results
save(envelope_MAT_file, 'outliers', '-append');