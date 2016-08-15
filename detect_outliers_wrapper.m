function detect_outliers(schalk_subjid, varargin)

% number of standard deviations above the mean to exclude
outlier_threshold = 6;

% general-purpose ecog analysis code
addpath([root_directory '/ecog-analysis-code']);

% directory to save results to
analysis_directory = [root_directory '/naturalsound-ecog/analysis/gamma-stimulus-response-from-preproc'];
if ~exist(analysis_directory, 'dir');
    mkdir(analysis_directory);
end

% load data
load([analysis_directory '/' schalk_subjid '.mat'], 'gamma_stimulus_response');
[n_smps, n_stimuli, n_runs, n_electrodes] = size(gamma_stimulus_response);
D = reshape(gamma_stimulus_response, [n_smps*n_stimuli, n_runs, n_electrodes]);

% determine cutoff for each electrode and run
mean_sd = quantile(D,[0.5 0.8413],1);
cutoffs = mean_sd(1,:,:) + diff(mean_sd) * outlier_threshold;

% use cutoff to detect outliers
outliers = D > repmat(cutoffs, [n_smps*n_stimuli, 1, 1]);
outliers = reshape(outliers, [n_smps, n_stimuli, n_runs, n_electrodes]);
save([analysis_directory '/' schalk_subjid '.mat'], 'outliers', '-append');
