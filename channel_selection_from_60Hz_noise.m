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
% 
% 2019-06-12 - Added functionality to compute 60 Hz noise within arrays

I.bw = 0.6;
I.frac = 0.2; % fraction of electrodes used to calculate standard deviation
I.thresh = 5; % threshold needed to calculate
I.array_inds = ones(1,size(signal,2)); % indices for different arrays
I.min_nchannels = 10; % minimum number of channels needed to calculate distribution
I.chnames = {};
I = parse_optInputs_keyvalue(varargin, I);

n_channels = size(signal,2);

noNaN_channels = find(all(~isnan(signal)));
signal_noNaN = signal(:, noNaN_channels);

% measure 60 Hz power in rms units
[b,a] = iirpeak(60/(sr/2), I.bw/(sr/2));
rms60Hz = mean(sqrt(filter(b, a, signal_noNaN).^2),1);

% infer mode
% [N,bin_centers] = hist(noise60Hz_rms, 1000);
% [~,xi] = max(N);
% mode_60Hz = bin_centers(xi);

unique_arrays = unique(I.array_inds);
rms60Hz_zscore = nan(1, n_channels);
for i = 1:length(unique_arrays)
    
    array_chan = find(I.array_inds(noNaN_channels)==unique_arrays(i));
    
    % fraction of channels to use to calculate z-score
    frac_corresponding_to_min = I.min_nchannels / length(array_chan);
    frac = min(max(frac_corresponding_to_min, I.frac), 1);
    
    % zsore using samples nearest the median
    if frac < 1
        rms60Hz_zcore_within_array = zscore_using_central_samples(rms60Hz(array_chan), frac);
        rms60Hz_zscore(noNaN_channels(array_chan)) = rms60Hz_zcore_within_array;
    end
end

% select good channels
bad_channels = (abs(rms60Hz_zscore) > I.thresh) | ~noNaN_channels;
good_channels = find(~bad_channels);

% plots
figure;
bar(1:n_channels,rms60Hz_zscore);
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

% save a list of the good channel names
if ~isempty(I.chnames)
    good_chnames = I.chnames(good_channels);
    fid = fopen([figure_fname '_good_channel_names.txt'],'w');
    fprintf(fid,'%s\n', good_chnames{:});
    fclose(fid);
    
    % save a list of the bad channel names
    bad_chnames = I.chnames(setdiff(1:n_channels, good_channels));
    fid = fopen([figure_fname '_bad_channel_names.txt'],'w');
    fprintf(fid,'%s\n', bad_chnames{:});
    fclose(fid);
end

% print out where figures are saved to
fprintf('Saving figures to:\n%s\n', figure_fname);
drawnow;