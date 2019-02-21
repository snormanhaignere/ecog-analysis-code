function [MAT_file, R] = electrode_envelope_reliability_across_runs(exp, subjid, preproc_idstring, varargin)

% Calculates the test-retest reliability of the stimulus driven envelope
% responses. Assumes the same stimuli are presented each run.
%
% 2016-1-31: Created by Sam NH
%
% 2016-10-22: Integrated with preprocessing scripts, made a bit more general
% 
% 2019-01-21: Changed parameter handling

%% Parameters

% window to include when calculating the correlation
I.win = [];

% whether or not to remove the average response timecourse across all stimuli
% the goal is to test the reliability of stimulus specific effects
I.demean = false;

% whether or not to time-average the response within the window
% before calculating the correlation
I.timeav = false;

% duration in seconds of chunk size for permutation test
I.permdur = 0.5;

% baseline period used to calculate PSC values
I.basewin = 0.5;

% number of permutations for significance testing
I.nperms = 1000;

% whether to compute splithalf 
% alternative is to correlate every run with every other run
I.splithalf = true;

% whether or not to spearman brown correct
I.spearbrown = false;

% whether or to overwrite the existing results
I.overwrite = false;

% whether or not to plot results
I.plot = true;

% optionally overwrite defaults
[I, ~, C_value] = parse_optInputs_keyvalue(varargin, I);

% string identifying this analysis
reliability_idstring = [struct2string(I, 'include_fields', {'nperms', 'splithalf', 'spearbrown'}) ...
    '_' struct2string(C_value, 'omit_field', {'overwrite', 'plot'})];
if isempty(reliability_idstring); reliability_idstring = 'default'; end
if reliability_idstring(end)=='_'; reliability_idstring = reliability_idstring(1:end-1); end
param_idstring = [preproc_idstring '/' reliability_idstring];

%% Directories, file names

global root_directory;

% directory for this project
project_directory = [root_directory '/' exp];

% director to save analyses to
analysis_directory = [project_directory '/analysis/envelope-reliability/' subjid '/' param_idstring];
if ~exist(analysis_directory,'dir'); mkdir(analysis_directory); end

% directory to save plots to
figure_directory = strrep(analysis_directory, 'analysis', 'figures');
if ~exist(figure_directory,'dir'); mkdir(figure_directory); end

MAT_file = [analysis_directory '/reliability.mat'];

%% Plot median reliability for each electrode and stimulus, with or without demeaning the timecourse first

if ~exist(MAT_file, 'file') || I.overwrite
    
    % load the envelope data
    % average reps within runs, concatenate across runs
    % D, outiers:
    % -> samples/time x stimuli x runs x electrodes
    [D, t] = load_envelopes(exp, subjid, preproc_idstring, ...
        'basewin', I.basewin, 'interpNaNs', true);
    [n_smps, n_stimuli, n_runs, n_electrodes] = size(D);
    env_sr = 1/(t(2)-t(1));
    
    % optionally remove average response across stimuli
    if I.demean
        D = bsxfun(@minus, D, nanmean( nanmean(D, 2), 3));
    end
    
    % optionally select subset of the response window
    if ~isempty(I.win)
        xi = t >= I.win(1) & t <=I.win(2);
        D = D(xi,:,:,:);
        n_smps = size(D,1);
        t = t(xi); %#ok<NASGU>
    end
    
    % optionally average across time
    if I.timeav
        D = squeeze_dims(nanmean(D,1),1);
    else
        D = reshape(D,[n_smps * n_stimuli, n_runs, n_electrodes]);
    end
    
    % compute reliability
    R = reliability(D, 'removeNaNs', true, ...
        'splithalf', I.splithalf, 'spearbrown', I.spearbrown, 'nperms', I.nperms, ...
        'chunksize', round(I.permdur*env_sr));
    
    % save results
    save(MAT_file, 'R');
    
else
    
    load(MAT_file, 'R');
    
end

%% Plot

if ~I.plot
    return;
end

if ~exist('R', 'var')
    load(MAT_file, 'R');
end

stat_names = {'corr', 'logP-gauss', 'logP-counts', 'zstat'};
for i = 1:length(stat_names)
    
    switch stat_names{i}
        case 'corr'
            stat_to_plot = R.corr;
            ylim_bounds = [0 1];
            ylabel_str = 'Test-Retest Correlation (r)';
        case 'logP-gauss'
            stat_to_plot = R.logP_gauss;
            ylim_bounds = [0 10];
            ylabel_str = 'Significance (-log10[p])';
        case 'logP-counts'
            stat_to_plot = R.logP_counts;
            ylim_bounds = [0 -log10(1/I.nperms)];
            ylabel_str = 'Significance (-log10[p])';
        case 'zstat'
            stat_to_plot = R.z;
            ylim_bounds = [0 10];
            ylabel_str = 'Significance (Z)';
        otherwise
            error('No matching case for stat_name %s\n', stat_names{i});
    end
    
    plot_electrode_statistic(stat_to_plot, stat_names{i}, 'ylim', ylim_bounds);
    ylabel(ylabel_str);
    
    export_fig([figure_directory '/' stat_names{i} '.pdf'],'-pdf','-transparent');
    export_fig([figure_directory '/' stat_names{i} '.png'],'-png','-transparent','-r100');
    
    %     % plot median reliability across stimuli
    %     n_electrodes_per_row = 50;
    %     total_num_rows = ceil(n_electrodes / n_electrodes_per_row);
    %     figure;
    %     set(gcf,'Position',[0 0 1200 1200]);
    %     subplot(total_num_rows,1,1);
    %     bar(1:n_electrodes,  stat_to_plot, 'FaceColor', [0.5 0.5 0.5]);
    %     xlim([0 n_electrodes+1]);
    %     ylabel(ylabel_str);
    %     xlabel('Research Electrode Number');
    %     for j = 1:ceil(n_electrodes / n_electrodes_per_row)
    %         subplot(total_num_rows+1,1,j+1);
    %
    %         electrodes_to_plot = (1:n_electrodes_per_row) + n_electrodes_per_row*(j-1);
    %         xi = electrodes_to_plot<=n_electrodes;
    %         electrodes_to_plot = electrodes_to_plot(xi);
    %
    %         bar(1:length(electrodes_to_plot),  stat_to_plot(electrodes_to_plot), 'FaceColor', [0.5 0.5 0.5]);
    %         xlim([0 n_electrodes_per_row+1]);
    %         ylim(ylim_bounds);
    %         grid on;
    %         set(gca, 'XTick', 1:length(electrodes_to_plot), 'XTickLabel',electrodes_to_plot);
    %     end

end

% % scatter plot comparison of responses with and without demeaning
% figure;
% plot(R.corr, reliability_after_demeaning, 'ko', 'LineWidth', 2);
% bounds = [min([R.corr,reliability_after_demeaning]),1];
% xlim(bounds); ylim(bounds);
% hold on; plot(bounds, bounds, 'r--', 'LineWidth',2);
% xlabel('Reliability (Without Demeaning)'); ylabel('Reliability (With Demeaning)');
% export_fig([figure_directory '/' 'timecourse-reliability-effect-of-demeaning-' num2str(n_runs) 'runs' '.pdf'],'-pdf','-transparent');
% export_fig([figure_directory '/' 'timecourse-reliability-effect-of-demeaning-' num2str(n_runs) 'runs' '.png'],'-png','-transparent','-r100');





