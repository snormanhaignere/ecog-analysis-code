function [D, si, ei, stim_names, n_runs, MAT_file] = format_ecog_data_standardized(exp, all_subjid, varargin)

% Formats data from all subjects into a matrix used by other programs
%
% The matrix D is formatted as time x stimuli x repetitions x electrodes
%
% Data from all subjects is combined, and subject indices for each electrode are
% stored as a vector (si). The corresponding ids are stored in all_subjid.

%% Directories

% directory with projects
global root_directory;

% outide code directories
% addpath([root_directory '/ecog-analysis-code']);
addpath([root_directory '/general-analysis-code']);
addpath([root_directory '/export_fig_v2']);

% experiment name and all subject IDs
% exp = 'naturalsound-ecog';
% all_subjid = {'AMC045', 'AMC047', 'AMC056', 'AMC062', 'AMC071', 'AMC078', 'AMC079', 'AMC081'};
n_subjects = length(all_subjid);

% preprocessing parameters
I.preproc_resp_win = [-1 3];
preprocessing_analysis_name = 'standard_preproc_params';
freq_cutoffs = [70 140]';
env_sr = 100;

% reliability threshold used to select data
I.reliability_threshold = NaN;

% final sampling rate
I.analysis_sr = 100;

% output time window
I.time_win = [];

% whether or not to average across runs
I.average_runs = false;

% whether or not to overwrite an existing MAT file
I.overwrite = false;

% whether or not to enter debug keyboard mode
I.keyboard = false;

I = parse_optInputs_keyvalue(varargin, I);

% default timewindow if not specified is set based on the preprocessing
% time window
if isempty(I.time_win)
    I.time_win = [0, I.preproc_resp_win(2) - 1/I.analysis_sr];
end

% enter debug mode
if I.keyboard
    keyboard;
end

preproc_idstring = [...
    preprocessing_analysis_name '_' ...
    num2str(freq_cutoffs(1)) '-' num2str(freq_cutoffs(2)) 'Hz' ...
    '_prewin' num2str(I.preproc_resp_win(1)) 'to' num2str(I.preproc_resp_win(2))];

% file to save results to
output_directory = [root_directory '/' exp '/analysis/formatted-data'];
if ~exist(output_directory, 'dir'); mkdir(output_directory); end
MAT_file = [output_directory '/'  sprintf('%s_', all_subjid{:}) ...
    preproc_idstring '_reliability' num2str(I.reliability_threshold) ...
    '_' num2str(I.analysis_sr) 'Hz' ...
    '_finalwin' num2str(I.time_win(1)) 'to' num2str(I.time_win(2)) ...
    '_avruns' num2str(I.average_runs) '_newsubjorder.mat'];

% load MAT file and return if already created
if exist(MAT_file, 'file') && ~I.overwrite
    load(MAT_file, 'D', 'si', 'ei', 'n_runs', 'stim_names');
    return;
end

%% Reliability information

if ~isnan(I.reliability_threshold)
    
    reliable_electrodes = cell(1,n_subjects);
    n_reliable_electrodes = nan(1,n_subjects);
    rel_corr = cell(1,n_subjects);
    for i = 1:n_subjects
        
        % load test-retest correlation
        reliability_MAT_file = electrode_envelope_reliability_across_runs(exp, all_subjid{i}, ...
            preprocessing_analysis_name, freq_cutoffs, I.preproc_resp_win, 'win', [0 3], 'plot', false);
        load(reliability_MAT_file,'reliability_corr');
        
        reliable_electrodes{i} = find(reliability_corr > I.reliability_threshold);
        rel_corr{i} = reliability_corr(reliability_corr > I.reliability_threshold);
        
        % remove errant electrode
        if strcmp(all_subjid{i}, 'AMC071')
            [~,xi] = setdiff(reliable_electrodes{i}, 138);
            reliable_electrodes{i} = reliable_electrodes{i}(xi);
            rel_corr{i} = rel_corr{i}(xi);
            clear xi;
        end
        if strcmp(all_subjid{i}, 'AMC062')
            [~,xi] = setdiff(reliable_electrodes{i}, 186);
            reliable_electrodes{i} = reliable_electrodes{i}(xi);
            rel_corr{i} = rel_corr{i}(xi);
            clear xi;
        end
        
        % store number of reliable electrodes
        n_reliable_electrodes(i) = length(reliable_electrodes{i});
        
    end
    
end

%% Envelope data

% load data
D_cell = cell(1, n_subjects);
electrode_indices = cell(1, n_subjects);
n_runs = nan(1, n_subjects);
n_electrodes = nan(1, n_subjects);
for i = 1:n_subjects
    
    [D_cell{i}, t, stim_names] = load_envelopes(exp, all_subjid{i}, ...
        preprocessing_analysis_name, freq_cutoffs, I.preproc_resp_win);
    assert(abs(1/(t(2)-t(1)) - env_sr)<1e-10);
    
    % select reliable electrodes
    if ~isnan(I.reliability_threshold)
        D_cell{i} = D_cell{i}(:,:,:,reliable_electrodes{i});
        electrode_indices{i} = reliable_electrodes{i};
    else
        electrode_indices{i} = 1:size(D_cell{i},4);
    end
    
    % dimensions
    [n_smps_per_stim, n_stimuli, n_runs(i), n_electrodes(i)] = size(D_cell{i});
    
end

%% Further preprocessing

% initialize
if I.average_runs
    D = nan(n_smps_per_stim, n_stimuli, sum(n_electrodes));
else
    D = nan(n_smps_per_stim, n_stimuli, max(n_runs), sum(n_electrodes));
end
si = nan(1, sum(n_electrodes));
ei = nan(1, sum(n_electrodes));

% average runs and unwrap
for i = 1:n_subjects
    xi = (1:n_electrodes(i)) + sum(n_electrodes(1:i-1));
    if I.average_runs
        D(:,:,xi) = squeeze_dims(nanmean(D_cell{i},3),3);
    else
        D(:,:,1:n_runs(i),xi) = D_cell{i};
    end
    si(xi) = i;
    ei(xi) = electrode_indices{i};
    clear xi;
end

% interpolate NaN values
D = interpNaN_ndarray(1:n_smps_per_stim, D, 'pchip');

% resample to desired window and sample rate
D = resample_and_window(D, I.preproc_resp_win, env_sr, I.time_win, I.analysis_sr);

%% Save results

save(MAT_file, 'D', 'ei', 'si', 'n_runs', 'stim_names', '-v7.3');
fprintf('Results saved here:\n%s\n', MAT_file); drawnow;

%% Compute the variance for each electrode

% electrode_std_cell = cell(size(D_cell));
% for i = 1:n_subjects
%
%     X = reshape(D_cell{i}, [n_smps_per_stim*n_stimuli, n_runs(i), n_electrodes(i)]);
%
%     %     % demean
%     %     X = bsxfun(@minus, X, nanmean(X,1));
%     %
%     %     % standardize
%     %     Z = bsxfun(@times, nanmean(nanstd(X,[],1),2), 1./nanstd(X,[],1));
%     %     X = bsxfun(@times, X, Z);
%     %     clear Z;
%
%     % all pairs of runs
%     pairs = nchoosek(1:n_runs(i),2);
%     n_pairs = size(pairs,1);
%
%     % estimated error for a single measurement based on all pairs of runs
%     E = nan(size(X,1),n_pairs,n_electrodes(i));
%     for j = 1:n_pairs
%         E(:,j,:) = X(:,pairs(j,1),:) - X(:,pairs(j,2),:);
%     end
%
%     % error for the average of all measurements
%     % sqrt(2) accounts for the fact that variances add
%     % sqrt(n_repeated_measurements) accounts for standard error calculation
%     electrode_std_cell{i} = squeeze(sqrt(nanmean(nanmean(E.^2,1),2))/sqrt(2)); %*sqrt(n_repeated_measurements));
%
% end

%% Further preprocessing
% 
% % average runs and unwrap
% D = [];
% si = [];
% for i = 1:n_subjects
%     D = cat(3, D, squeeze_dims(nanmean(D_cell{i},3),3));
%     si = [si, i*ones(1, size(D_cell{i},4))];
% end
% electrode_std = cat(1, electrode_std_cell{:})';
% 
% % interpolate NaN values
% D = interpNaN_ndarray(1:n_smps_per_stim, D, 'pchip');
% 
% % resample to desired window and sample rate
% D = resample_and_window(D, I.preproc_resp_win, env_sr, [0 2.96], 25);
% 
% delete([root_directory '/naturalsound-ecog/analysis/formatted-data/reliability0.1.mat']);
% save([root_directory '/naturalsound-ecog/analysis/formatted-data/reliability0.1.mat'], 'D', 'si', 'electrode_std', '-v7.3');
% 
% 
% ei = cat(2, reliable_electrodes{:})';
% save([root_directory '/naturalsound-ecog/analysis/formatted-data/reliability0.1.mat'], '-append', 'ei')
% save([root_directory '/naturalsound-ecog/analysis/formatted-data/reliability0.1.mat'], '-append', 'all_subjid')
% 
% %%
% 
% 
% % load data
% D_cell = cell(1, n_subjects);
% n_electrodes = nan(1, n_subjects);
% for i = 1:n_subjects
%     
%     [D_cell{i},t] = load_envelopes(exp, all_subjid{i}, ...
%         preprocessing_analysis_name, freq_cutoffs, I.preproc_resp_win);
%     
%     % dimensions
%     [n_smps_per_stim, n_stimuli, n_runs(i), n_electrodes(i)] = size(D_cell{i});
%     
% end
% 
% % average runs and unwrap
% D = nan(n_smps_per_stim, n_stimuli, 7, sum(n_electrodes));
% si = [];
% for i = 1:n_subjects
%     dims = size(D_cell{i});
%     xi = (1:n_electrodes(i)) + sum(n_electrodes(1:i-1));
%     D(:,:,1:dims(3),xi) = D_cell{i}(:, :, :, :);
%     si = [si, i*ones(1, size(D_cell{i},4))];
% end
% 
% % interpolate NaN values
% D = interpNaN_ndarray(1:n_smps_per_stim, D, 'pchip');
% 
% % resample to desired window and sample rate
% % Y = resample_and_window(X, I.preproc_resp_win, 100, [0 2.99], 100);
% 
% D = D(101:400,:,:,:);
% 
% 
% % delete([root_directory '/naturalsound-ecog/analysis/formatted-data/all-electrodes-and-runs-100Hz-0-3.mat']);
% save([root_directory '/naturalsound-ecog/analysis/formatted-data/all-electrodes-and-runs-100Hz-0-3.mat'], 'D', 'si', '-v7.3');
% 
% 
% %% Create splits
% 
% % average runs and unwrap
% D1 = [];
% for i = 1:n_subjects
%     D1 = cat(3, D1, squeeze_dims(nanmean(D_cell{i}(:,:,1:2:end,:),3),3));
% end
% 
% % interpolate NaN values
% D1 = interpNaN_ndarray(1:n_smps_per_stim, D1, 'pchip');
% 
% % resample to desired window and sample rate
% D1 = resample_and_window(D1, I.preproc_resp_win, 100, [0 2.96], 25);
% 
% % average runs and unwrap
% D2 = [];
% for i = 1:n_subjects
%     D2 = cat(3, D2, squeeze_dims(nanmean(D_cell{i}(:,:,2:2:end,:),3),3));
% end
% electrode_std = cat(1, electrode_std_cell{:})';
% 
% % interpolate NaN values
% D2 = interpNaN_ndarray(1:n_smps_per_stim, D2, 'pchip');
% 
% % resample to desired window and sample rate
% D2 = resample_and_window(D2, I.preproc_resp_win, 100, [0 2.96], 25);
% 
% delete([root_directory '/naturalsound-ecog/analysis/formatted-data/reliability0.1-splits.mat']);
% save([root_directory '/naturalsound-ecog/analysis/formatted-data/reliability0.1-splits.mat'], 'D1', 'D2', '-v7.3');
% 
% 
% xi = 20;
% A = D(:,:,xi);
% bounds = quantile(A(:), 0.99);
% subplot(1,3,1);
% imagesc(D(:,:,xi), [-bounds bounds]);
% subplot(1,3,2);
% imagesc(D1(:,:,xi), [-bounds bounds]);
% subplot(1,3,3);
% imagesc(D2(:,:,xi), [-bounds bounds]);
% 
% %% Training and testing stimuli
% 
% load([root_directory '/' exp '/code/category_regressors.mat'],'C');
% 
% rand_seed = 0;
% ResetRandStream2(0);
% partition_fractions = [0.6,0.2,0.2];
% n_partition_fractions = length(partition_fractions);
% partition_index = nan(length(C.category_assignments),1);
% for i = 1:length(C.category_labels);
%     
%     xi = find(C.category_assignments == i);
%     xi = xi(randperm(length(xi)));
%     
%     count = 0;
%     for j = 1:length(partition_fractions)
%         
%         start_ind = 1+count;
%         if j == length(partition_fractions)
%             end_ind = length(xi);
%         else
%             end_ind = round(partition_fractions(j)*length(xi))+count;
%             end_ind = min(end_ind, length(xi));
%         end
%         yi = start_ind:end_ind;
%         
%         partition_index(xi(yi)) = j-1;
%         count = count + length(yi);
%     end
% end
% 
% save([root_directory '/naturalsound-ecog/analysis/formatted-data/partition' ...
%     sprintf('-%.1f', partition_fractions) '-randseed' num2str(rand_seed) '.mat'],...
%     'partition_index', '-v7.3');


