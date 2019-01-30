function [good_channels, rms60Hz] = ...
    channel_selection_from_60Hz_noise(signal, sr, figure_fname, varargin)

% Uses 60 Hz noise to detect good/bad channels
% 
% Derived from Schalk lab code
% 
% signal: [time x electrode] singal matrix
% 
% sr: sampling rate of the signal
% 
% figure_fname: name of file to save plots to
% 
% 2016-08-12 - Created, Sam NH
% 
% 2016-09-20 - Modified to use the mode instead of the median to detect outliers
% 
% 2016-09-23 - Switched back to median, but use alternative z-scoring procedure
% to detect bad electrodes, Sam NH

I.bw = 0.6;
I.frac = 0.2; % fraction of electrodes used to calculate standard deviation
I.thresh = 5; % threshold needed to calculate
I = parse_optInputs_keyvalue(varargin, I);

n_channels = size(signal,2);

noNaN_channels = all(~isnan(signal));
signal_noNaN = signal(:, noNaN_channels);

% measure 60 Hz power in rms units
[b,a] = iirpeak(60/(sr/2), I.bw/(sr/2));
rms60Hz = mean(sqrt(filter(b, a, signal_noNaN).^2),1);

% infer mode
% [N,bin_centers] = hist(noise60Hz_rms, 1000);
% [~,xi] = max(N);
% mode_60Hz = bin_centers(xi);

% zsore using samples nearest the median
rms60Hz_zcore = zscore_using_central_samples(rms60Hz, I.frac);
rms60Hz_zcore = fillin_NaN(rms60Hz_zcore, noNaN_channels, 2);

% select good channels
good_channels = find(abs(rms60Hz_zcore) < I.thresh);

% plots
figure;
bar(1:n_channels,rms60Hz_zcore);
hold on;
plot([1,n_channels], I.thresh*[1 1], 'r--');
plot([1,n_channels], -I.thresh*[1 1], 'r--');
xlabel('Electrodes'); ylabel('RMS Power (Dev. from Median)')
box off;
export_fig([figure_fname '.pdf'], '-pdf', '-transparent');
close all;

% save a list of the good channels
fid = fopen([figure_fname '_good_channels.txt'],'w');
fprintf(fid,'%d\n', good_channels);
fclose(fid);

% save a list of the bad channels
fid = fopen([figure_fname '_bad_channels.txt'],'w');
fprintf(fid,'%d\n', setdiff(1:n_channels, good_channels));
fclose(fid);