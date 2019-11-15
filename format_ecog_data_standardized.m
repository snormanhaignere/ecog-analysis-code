function [D, t, si, ei, stim_names, n_runs, relstat, MAT_file, param_idstring] = format_ecog_data_standardized(...
    exp, all_subjid, preproc_idstring, varargin)

% Formats data from all subjects into a matrix used by other programs
%
% The matrix D is formatted as time x stimuli x repetitions x electrodes
%
% Data from all subjects is combined, and subject indices for each electrode are
% stored as a vector (si). The corresponding ids are stored in all_subjid.
% 
% 2019-01-21: Changed parameter handling, Sam NH

%% Directories

% directory with projects
global root_directory;

% outide code directories
% addpath([root_directory '/ecog-analysis-code']);
addpath([root_directory '/general-analysis-code']);
addpath([root_directory '/export_fig_v3']);

% experiment name and all subject IDs
n_subjects = length(all_subjid);

% reliability threshold used to select data
I.relstat = 'corr';
I.relthresh = NaN;
I.splithalf = true;
I.spearbrown = false;
I.permdur = 0.5;

% final sampling rate
I.sr = 100;

% output time window
I.win = [];

% baseline period used to calculate PSC values
I.basewin = 0.5;

% whether or not to average across runs
I.avruns = false;

% whether or not to overwrite an existing MAT file
I.overwrite = false;

% whether or not to interpolate NaN values
I.interpNaNs = true;

% whether or not to hand exclude certain electrodes/runs
I.handexclude = false;

% whether or not to enter debug keyboard mode
I.keyboard = false;

[I, ~, C_value] = parse_optInputs_keyvalue(varargin, I);

% enter debug mode
if I.keyboard
    keyboard;
end

% idstring to identify this analysis
format_idstring = struct2string(C_value, 'omit_field', {'overwrite', 'keyboard'});
if isempty(format_idstring); format_idstring = 'default-format'; end
subj_string = sprintf('%s_', all_subjid{:}); subj_string = subj_string(1:end-1);
param_idstring = [preproc_idstring '/' subj_string '/' format_idstring];
clear subj_string;

% file to save results to
output_directory = [root_directory '/' exp '/analysis/formatted-data/' param_idstring];
if ~exist(output_directory, 'dir'); mkdir(output_directory); end
MAT_file = [output_directory '/data.mat'];

% load MAT file and return if already created
if exist(MAT_file, 'file') && ~I.overwrite
    load(MAT_file, 'D', 't', 'si', 'ei', 'n_runs', 'stim_names', 'relstat');
    return;
end

%% Reliability information

if ~isnan(I.relthresh)
    
    reliable_electrodes = cell(1,n_subjects);
    n_reliable_electrodes = nan(1,n_subjects);
    relstat_cell = cell(1,n_subjects);
    for i = 1:n_subjects

        % load test-retest correlation
        [~, R] = electrode_envelope_reliability_across_runs(....
            exp, all_subjid{i}, preproc_idstring, 'win', I.win, 'plot', false, ...
            'splithalf', I.splithalf, 'spearbrown', I.spearbrown, 'permdur', I.permdur);
        
        reliable_electrodes{i} = find(R.(I.relstat) > I.relthresh);
        relstat_cell{i} = R.(I.relstat)(R.(I.relstat) > I.relthresh);
        
        % remove errant electrode
        if strcmp(all_subjid{i}, 'AMC071')
            [~,xi] = setdiff(reliable_electrodes{i}, 138);
            reliable_electrodes{i} = reliable_electrodes{i}(xi);
            relstat_cell{i} = relstat_cell{i}(xi);
            clear xi;
        end
        if strcmp(all_subjid{i}, 'AMC062')
            [~,xi] = setdiff(reliable_electrodes{i}, 186);
            reliable_electrodes{i} = reliable_electrodes{i}(xi);
            relstat_cell{i} = relstat_cell{i}(xi);
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
    
    [D_cell{i}, t, stim_names, electrode_research_numbers] = ...
        load_envelopes(exp, all_subjid{i}, preproc_idstring, ...
        'basewin', I.basewin, 'interpNaNs', I.interpNaNs, ...
        'win', I.win, 'sr', I.sr);
    
    % select reliable electrodes
    if ~isnan(I.relthresh)
        D_cell{i} = D_cell{i}(:,:,:,reliable_electrodes{i});
        electrode_indices{i} = electrode_research_numbers(reliable_electrodes{i});
    else
        electrode_indices{i} = electrode_research_numbers;
    end
    
    % dimensions
    [n_smps_per_stim, n_stimuli, n_runs(i), n_electrodes(i)] = size(D_cell{i});
    
end

%% Exclude electrodes

if I.handexclude
    load([root_directory '/' exp '/analysis/hand-exclude.mat'], 'exclude');
    for i = 1:size(exclude,1)
        xi = ismember(all_subjid, exclude{i,1});
        assert(sum(xi)<2);
        if sum(xi)>0
            yi = electrode_indices{xi}==exclude{i,2};
            assert(sum(yi)==1);
            if length(exclude{i,3})==1 && isnan(exclude{i,3})
                fprintf('Excluding electrode %d\n', exclude{i,2});
                D_cell{xi} = D_cell{xi}(:,:,:,~yi);
                electrode_indices{xi} = electrode_indices{xi}(~yi);
                n_electrodes(xi) = n_electrodes(xi)-1;
            else
                D_cell{xi}(:,:,exclude{i,3},yi) = NaN;
            end
        end
        clear xi yi;
    end
end

%% Unwrap to matrix

% initialize
if I.avruns
    D = nan(n_smps_per_stim, n_stimuli, sum(n_electrodes));
else
    D = nan(n_smps_per_stim, n_stimuli, max(n_runs), sum(n_electrodes));
end
si = nan(1, sum(n_electrodes));
ei = nan(1, sum(n_electrodes));
if ~isnan(I.relthresh)
    relstat = nan(1, sum(n_electrodes));
else
    relstat = [];
end

% unwrap across subjects while storing electrode and subject indices
% optioanlly average across runs
for i = 1:n_subjects
    xi = (1:n_electrodes(i)) + sum(n_electrodes(1:i-1));
    if I.avruns
        D(:,:,xi) = squeeze_dims(nanmean(D_cell{i},3),3);
    else
        D(:,:,1:n_runs(i),xi) = D_cell{i};
    end
    si(xi) = i;
    ei(xi) = electrode_indices{i};
    if ~isnan(I.relthresh)
        relstat(xi) = relstat_cell{i};
    end
    clear xi;
end

%% Save results

save(MAT_file, 'D', 't', 'ei', 'si', 'n_runs', 'stim_names', 'relstat', '-v7.3');
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
% D = resample_and_window(D, preproc_resp_win, I.preproc_sr, [0 2.96], 25);
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
%         preprocessing_analysis_name, freq_cutoffs, preproc_resp_win);
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
% % Y = resample_and_window(X, preproc_resp_win, 100, [0 2.99], 100);
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
% D1 = resample_and_window(D1, preproc_resp_win, 100, [0 2.96], 25);
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
% D2 = resample_and_window(D2, preproc_resp_win, 100, [0 2.96], 25);
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


%%



