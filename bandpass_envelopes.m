function env = bandpass_envelopes(signal, signal_sr, env_sr, ...
    band_in_Hz, filter_order, figure_directory, figure_fname_prefix, ...
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
% band_in_Hz: 2-dimensional vector of frequency cutoffs
% 
% filter_orders: order of the butterworth filters
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
    
% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', filter_order, ...
    band_in_Hz(1), band_in_Hz(2), signal_sr);
Hd = design(h,'butter');
[B, A] = sos2tf(Hd.sosMatrix,Hd.scaleValues);

% apply filter
fprintf('Applying filter...\n');
subb = filtfilt(B,A,signal);

% measure envelope
fprintf('Calculating envelopes...\n');
env = abs(hilbert(subb));

% plots
fname = [figure_directory '/' figure_fname_prefix];
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

