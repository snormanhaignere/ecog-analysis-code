function signal = hp_filt(signal, P, sr)

% High-pass filter ECoG signal to remove very low frequencies (e.g. < 0.5 Hz)
% 
% Derived from Schalk lab code
% 
% signal: [time x electrode] singal matrix
% 
% P: parameters, see preprocessing_parameters.m
% 
% sr: sampling rate of the signal
% 
% 2016-08-12 - Created, Sam NH

% design filter
Hd = design(...
    fdesign.highpass('N,F3dB', P.hp_filt_order, P.hp_filt_cutoff_in_Hz, sr), ...
    'butter');

% convert to pole form
[b, a] = sos2tf(Hd.sosMatrix,Hd.scaleValues);

% apply filter
signal = filtfilt(b,a,signal);