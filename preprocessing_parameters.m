function P = preprocessing_parameters

% bandwidth of the peak filter used to measure 60 Hz noise
% specifies the 3 dB down point in Hz
% parameter value derived from schalk lab code
P.bw_60Hz_peak_filt = 0.6;

% filter order and cutoff of the highpass filter
% used to remove very low frequencies
% parameters values derived from schalk lab code
P.hp_filt_order = 4;
P.hp_filt_cutoff_in_Hz = 0.5;

% parameters of the notch filter used to remove 60 Hz noise and its harmonics
P.notch_n_harmonics = 4;
P.notch_bw = 1; % 3 dB down/up in Hz

% parameter of bandpass filters
P.bandpass_env_sr = 100; 
P.bandpass_cutoffs_in_Hz = [70 140]'; % cutoffs
P.bandpass_filter_orders = 6 * ones(1,3); % order of filters, one per bandwidth

% threshold used to detect outliers in envelopes
% the threshold is sort of in units of standard deviation
% but the standard deviation is measured using the difference between the median
% of the distribution and the 84th percential, which for a gaussian would give a
% single standard deviation, the advantage of this measure is that it's less
% sensitive to outliers
P.envelope_outlier_threshold = 6; 