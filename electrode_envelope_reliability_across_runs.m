function MAT_file = electrode_envelope_reliability_across_runs(exp, subjid, ...
    analysis_name, freq_cutoffs, resp_win, varargin)

% Calculates the test-retest reliability of the stimulus driven envelope
% responses. Assumes the same stimuli are presented each run.
%
% 2016-1-31: Created by Sam NH
%
% 2016-10-22: Integrated with preprocessing scripts, made a bit more general

%% Parameters

% window to include when calculating the correlation
I.win = resp_win;

% whether or not to remove the average response timecourse across all stimuli
% the goal is to test the reliability of stimulus specific effects
I.demean = false;

% whether or not to time-average the response within the window
% before calculating the correlation
I.time_average = false;

% whether or to overwrite the existing results
I.overwrite = false;

% duration in seconds of chunk size for permutation test
I.perm_chunk_duration = 0.5;

% whether or not to plot results
I.plot = true;

% optionally overwrite defaults
I = parse_optInputs_keyvalue(varargin, I);

%% Directories, file names

global root_directory;

% directory for this project
project_directory = [root_directory '/' exp];

% director to save analyses to
analysis_directory = [project_directory ...
    '/analysis/electrode-envelope-reliability/' subjid];
if ~exist(analysis_directory,'dir')
    mkdir(analysis_directory);
end

% directory to save plots to
figure_directory = [project_directory ...
    '/figures/electrode-envelope-reliability/' subjid];
if ~exist(figure_directory,'dir')
    mkdir(figure_directory);
end

% string identifying the frequency range
freq_string = [num2str(freq_cutoffs(1)) '-' num2str(freq_cutoffs(2)) 'Hz'];

% mat files to save results to
idstring = [freq_string '_' analysis_name ...
    '_win' num2str(I.win(1)) '_' num2str(I.win(2))];

% additional flags
if I.time_average
    idstring = [idstring '_time_averaged'];
end
if I.demean
    idstring = [idstring '_demeaned'];
end

MAT_file = [analysis_directory '/' idstring '.mat'];

%% Plot median reliability for each electrode and stimulus, with or without demeaning the timecourse first


if ~exist(MAT_file, 'file') || I.overwrite
    
    % load the envelope data
    % average reps within runs, concatenate across runs
    % D, outiers:
    % -> samples/time x stimuli x runs x electrodes
    runs = read_runs(exp, subjid);
    n_runs = length(runs);
    for i = 1:n_runs
        
        input_MAT_file = [project_directory '/analysis/envelopes/' ...
            subjid '/r' num2str(runs(i)) '/envelopes_' freq_string '_' ...
            analysis_name '_mapped2stim_win' ...
            num2str(resp_win(1)) 'to' num2str(resp_win(2)) '.mat'];
        load(input_MAT_file, 'envelopes_mapped_to_stim', 'env_sr', ...
            'outliers_mapped_to_stim', 'stim_names');
        
        % average across any repetitions with a run
        if i == 1
            [n_smps, n_stimuli, ~, n_electrodes] = size(envelopes_mapped_to_stim);
            D = nan([n_smps, n_stimuli, n_runs, n_electrodes]);
        end
        try
            outliers_mapped_to_stim(isnan(outliers_mapped_to_stim)) = 0;
            envelopes_mapped_to_stim(logical(outliers_mapped_to_stim)) = NaN; %#ok<AGROW>
            D(:,:,i,:) = nanmean(envelopes_mapped_to_stim,3);
        catch
            keyboard
        end
        
    end
    
    % time-vector
    t = (0:n_smps-1)/env_sr + resp_win(1);
    
    % optionally remove average response across stimuli
    if I.demean
        D = bsxfun(@minus, D, nanmean( nanmean(D, 2), 3));
    end
    
    % optionally select subset of the response window
    if ~isequal(resp_win, I.win);
        xi = t >= I.win(1) & t <=I.win(2);
        D = D(xi,:,:,:);
        n_smps = size(D,1);
        t = t(xi); %#ok<NASGU>
    end
    
    % optionally average across time
    if I.time_average
        D = squeeze_dims(nanmean(D,1),1);
    else
        D = reshape(D,[n_smps * n_stimuli, n_runs, n_electrodes]);
    end
        
    % reliability per electrode
    reliability_corr = nan(1, n_electrodes);
    for i = 1:n_electrodes
        X = squeeze(D(:,:,i)); % samples/time x run for single electrode
        not_outliers = all(~isnan(X),2);
        r = corr(X(not_outliers,:)); % correlation of response timecourse for all pairs of runs
        average_r = mean(r(~eye(n_runs))); % average off diagonal entries
        reliability_corr(i) = correlation_power_analysis(average_r, n_runs, 0); % estimate reliability of entire dataset
    end
    
    % permutation test
    n_perms = 100;
    MAT_file_with_perms = strrep(MAT_file, '.mat', ['_perms' num2str(n_perms) '.mat']);
    if ~exist(MAT_file_with_perms,'file') || I.overwrite
        
        reliability_corr_null = nan(n_perms, n_electrodes);
        D_perm = D;
        for reliability_z = 1:n_perms
            
            if mod(reliability_z,10)==0
                fprintf('perm %d\n', reliability_z); drawnow;
            end
            
            % shuffle one-second chunks
            for j = 1:n_runs
                if I.time_average
                    D_perm(:,j,:) = D(randperm(n_stimuli),j,:);
                else
                    % divide into into chunks
                    chunk_size_smps = round(I.perm_chunk_duration * env_sr);
                    n_chunks = ceil(size(D,1)/chunk_size_smps);
                    xi = reshape(1:(chunk_size_smps * n_chunks), chunk_size_smps, n_chunks);
                    xi(size(D,1)+1:end) = NaN;
                    
                    % permute chunks
                    xi = xi(:,randperm(n_chunks));
                    xi = xi(:);
                    xi(isnan(xi)) = [];
                    
                    % assign based on permuted chunks
                    D_perm(:,j,:) = D(xi,j,:);
                end
            end
            
            % reliability per electrode
            for i = 1:n_electrodes
                X = squeeze(D_perm(:,:,i)); % samples/time x run for single stimulus
                not_outliers = all(~isnan(X),2);
                r = corr(X(not_outliers,:)); % correlation of response timecourse for all pairs of runs
                average_r = mean(r(~eye(n_runs))); % average off diagonal entries
                reliability_corr_null(reliability_z,i) = correlation_power_analysis(average_r, n_runs, 0); % estimate reliability of entire dataset
            end
            
        end
        
        save(MAT_file_with_perms, 'reliability_corr_null');
        
    else
        load(MAT_file_with_perms, 'reliability_corr_null');
    end
    
    % calculate significance
    m = mean(reliability_corr_null);
    sd = std(reliability_corr_null);
    reliability_z = (reliability_corr - m)./sd;
    reliability_logP = sign(reliability_z) .* -log10(2*normcdf(-abs(reliability_z), 0, 1));
    reliability_logP(isinf(reliability_logP)) = max(reliability_logP(~isinf(reliability_logP)));
    
    save(MAT_file,'reliability_corr', 'reliability_logP', 'reliability_z');
    
end

%% Plot

if ~I.plot
    return;
end

if ~exist('reliability_corr', 'var')
    load(MAT_file, 'reliability_corr');
end
if ~exist('reliability_logP', 'var')
    load(MAT_file, 'reliability_logP');
end
if ~exist('reliability_z', 'var')
    load(MAT_file, 'reliability_z');
end

stat_names = {'corr', 'logP', 'zstat'};
for i = 1:length(stat_names)
    
    switch stat_names{i}
        case 'corr'
            stat_to_plot = reliability_corr;
            ylim_bounds = [0 1];
            ylabel_str = 'Test-Retest Correlation (r)';
        case 'logP'
            stat_to_plot = reliability_logP;
            ylim_bounds = [0 10];
            ylabel_str = 'Significance (-log10[p])';
        case 'zstat'
            stat_to_plot = reliability_z;
            ylim_bounds = [0 10];
            ylabel_str = 'Significance (Z)';
        otherwise
            error('No matching case for stat_name %s\n', stat_names{i});
    end
    
    plot_electrode_statistic(stat_to_plot, stat_names{i}, 'ylim', ylim_bounds);
    ylabel(ylabel_str);
    
    export_fig([figure_directory '/' idstring '_' stat_names{i} '.pdf'],'-pdf','-transparent');
    export_fig([figure_directory '/' idstring '_' stat_names{i} '.png'],'-png','-transparent','-r100');
    
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
% plot(reliability_corr, reliability_after_demeaning, 'ko', 'LineWidth', 2);
% bounds = [min([reliability_corr,reliability_after_demeaning]),1];
% xlim(bounds); ylim(bounds);
% hold on; plot(bounds, bounds, 'r--', 'LineWidth',2);
% xlabel('Reliability (Without Demeaning)'); ylabel('Reliability (With Demeaning)');
% export_fig([figure_directory '/' 'timecourse-reliability-effect-of-demeaning-' num2str(n_runs) 'runs' '.pdf'],'-pdf','-transparent');
% export_fig([figure_directory '/' 'timecourse-reliability-effect-of-demeaning-' num2str(n_runs) 'runs' '.png'],'-png','-transparent','-r100');





