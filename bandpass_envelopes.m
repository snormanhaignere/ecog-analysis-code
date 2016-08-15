function envelopes = bandpass_envelopes(signal, signal_sr, env_sr, ...
    bands_in_Hz, filter_orders, figure_directory, figure_fname_prefix, ...
    good_channels, varargin)

% Measure envelopes from bandpass filters.
% 
% -- Inputs --
% 
% signal: [time x electrode] signal matrix
%
% signal_sr: sampling rate of the signal
% 
% env_sr: sampling rate to use for the envelopes (typically much lower than the
% signal sampling rate, e.g. 100 Hz vs. 1200 Hz)
% 
% bands_in_Hz: [2 x num_bands] matrix of cutoffs for the bandpass filters in Hz
% 
% filter_orders: order of the butterworth filters (default is 6), can specify a
% single filter order, in which case the filter order is applied to all filters,
% or a vector with length equal to the number of filters
% 
% figure_directory: directory to save figures to
% 
% figure_fname_prefix: file name prefix appended to the beginning of all saved
% files
% 
% good_channels: vector of channels deemed good, used for the purposes of
% plotting.
% 
% 
% -- Outputs --
% 
% envelopes: [time x bands x electrode] matrix with envelopes
% 
% 2016-1-26: Created by Sam NH

%% Setup

% default filter order if not specified
if nargin < 5
    filter_orders = 6;
end

% dimensionality of the signal
[n_signal_smps, n_electrodes] = size(signal);

% number of bands to extract
n_bands = size(bands_in_Hz,2);

% number of envelope samples
n_env_smps = ceil(n_signal_smps * env_sr / signal_sr);

if length(filter_orders) == 1
    filter_orders = filter_orders * ones(1,n_bands);
end

envelopes = nan(n_env_smps, n_bands, n_electrodes);
for i = 1:n_bands
    
    % filter frequencies outside of the gamma power range
    fprintf('Filtering for band %.2f - %.2f...\n', ...
        bands_in_Hz(1,i), bands_in_Hz(2,i)); drawnow;

    % Construct an FDESIGN object and call its BUTTER method.
    h  = fdesign.bandpass('N,F3dB1,F3dB2', filter_orders(i), bands_in_Hz(1,i), bands_in_Hz(2,i), signal_sr);
    Hd = design(h,'butter');
    [B, A] = sos2tf(Hd.sosMatrix,Hd.scaleValues);
    
    % apply filter
    fprintf('Applying filter...\n');
    subb = filtfilt(B,A,signal);
    
    % measure envelope
    fprintf('Calculating envelopes...\n');
    env = abs(hilbert(subb));
    
    % plot
    fname = [figure_directory '/' figure_fname_prefix '_bpfilt_' ...
        num2str(bands_in_Hz(1,i)) '-' num2str(bands_in_Hz(2,i)) 'Hz'];
    plot_example_spectra(signal, signal_sr, ...
        good_channels, [fname '_spectra_before']);
    plot_example_spectra(subb, signal_sr, ...
        good_channels, [fname '_spectra_after']);
    plot_example_timecourses(signal, signal_sr, ...
        good_channels, [fname '_timecourse_before']);
    plot_example_timecourses(subb, signal_sr, ...
        good_channels, [fname '_timecourse_after']);
    plot_example_timecourses(env, signal_sr, ...
        good_channels, [fname '_envelope']);
        
    % resample to desired rate
    fprintf('Downsampling...\n');
    env = resample(env, env_sr, signal_sr);
    
    % truncate
    env(env < 0) = 0;
    
    % store in matrix
    envelopes(:,i,:) = env;
    
end