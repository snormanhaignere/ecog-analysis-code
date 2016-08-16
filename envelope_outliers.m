function outliers = envelope_outliers(envelopes, outlier_threshold, figure_directory, figure_fname_prefix)

% Detects outliers in signal envelopes. The difference between the median and
% the 84th percential of the envelope distribution for each electrode is
% measured. For a Gaussian distribution this interval is equal to a standard
% deviation. Envelopes are considered outliers if they fall a certain number of
% 'standard deviations' above the median. 
% 
% 2016-08-15 - Created, Sam NH

% median and 84th percentile
% -> 2 x electrodes
median_sd = quantile(envelopes,[0.5 0.8413],1);

% cutoff
% -> 1 x electrodes
cutoffs = median_sd(1,:,:) + diff(median_sd) * outlier_threshold;

% detect outliers
% -> samples x electrodes
outliers = envelopes > repmat(cutoffs, [size(envelopes,1), 1]);

% plot
figure;
plot_electrode_statistic(mean(outliers)*100, '% Outliers');
box off;
fig_file = [figure_directory '/' figure_fname_prefix '_percent_outliers'];
set(gcf, 'PaperSize', [12 6]);
set(gcf, 'PaperPosition', [0.25 0.25 11.5 5.5]);
print([fig_file '.pdf'],'-dpdf');
print([fig_file '.png'],'-dpng', '-r100');
close all;