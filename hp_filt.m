function signal = hp_filt(signal, P, sr, varargin)

% High-pass filter ECoG signal to remove very low frequencies (e.g. < 0.5 Hz)
% 
% Derived from Schalk lab code
% 
% signal: [time x electrode] singal matrix
% 
% sr: sampling rate of the signal
% 
% 2016-08-12: Created, Sam NH
% 
% 2018-01-21: Changed parameter handling

I.order = 4;
I.cutoff = 0.5;
I = parse_optInputs_keyvalue(varargin, I);

% design filter
Hd = design(fdesign.highpass('N,F3dB', I.order, P.cutoff, sr), 'butter');

% convert to pole form
[b, a] = sos2tf(Hd.sosMatrix,Hd.scaleValues);

% apply filter
signal = filtfilt(b,a,signal);