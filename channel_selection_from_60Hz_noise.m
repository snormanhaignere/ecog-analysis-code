function [good_channels, noise60Hz_rms] = ...
    channel_selection_from_60Hz_noise(signal, P, sr, figure_fname)

% Uses 60 Hz noise to detect good/bad channels
% 
% Derived from Schalk lab code
% 
% signal: [time x electrode] singal matrix
% 
% P: parameters, see preprocessing_parameters.m
% 
% sr: sampling rate of the signal
% 
% figure_fname: name of file to save plots to
% 
% 2016-08-12 - Created, Sam NH
% 
% 2016-09-20 - Modified to use the mode instead of the median to detect outliers

n_channels = size(signal,2);

% measure 60 Hz power in rms units
[b,a] = iirpeak(60/(sr/2), P.bw_60Hz_peak_filt/(sr/2));
noise60Hz_rms = mean(sqrt(filter(b, a, signal).^2),1);

% infer mode
[N,bin_centers] = hist(noise60Hz_rms, 1000);
[~,xi] = max(N);
mode_60Hz = bin_centers(xi);

% measure the deviation from the mode
noise60Hz_deviation = noise60Hz_rms - mode_60Hz;

% threshold
dev_10percent = quantile(abs(noise60Hz_deviation),0.10);
dev_std = dev_10percent / (norminv(0.55,0,1) - norminv(0.45,0,1));
noise60Hz_thresh = 5*dev_std;
good_channels = find(abs(noise60Hz_deviation) < noise60Hz_thresh);

% plots
figure;
bar(1:n_channels,noise60Hz_deviation);
hold on;
plot([1,n_channels], noise60Hz_thresh*[1 1], 'r--');
plot([1,n_channels], -noise60Hz_thresh*[1 1], 'r--');
xlabel('Electrodes'); ylabel('RMS Power (Dev. from Median)')
box off;
export_fig([figure_fname '.pdf'], '-pdf', '-transparent');
close all;

% save a list of the good channels
fid = fopen([figure_fname '_good_channels.txt'],'w');
fprintf(fid,'%d\n', good_channels);
fclose(fid);