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

n_channels = size(signal,2);

% measure 60 Hz power in rms units
[b,a] = iirpeak(60/(sr/2), P.bw_60Hz_peak_filt/(sr/2));
noise60Hz_rms = mean(sqrt(filter(b, a, signal).^2),1);

% measure the deviation from the median
noise60Hz_deviation = noise60Hz_rms - median(noise60Hz_rms);

% threshold
noise60Hz_thresh = 10*mad(noise60Hz_deviation, 1);
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