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
P.bandpass_cutoffs_in_Hz = [5 10; 25 50; 70 140]'; % cutoffs
P.bandpass_filter_orders = 6; % order of butterworth filters
