function dh2 = EcogExtractFreqBand (d,infs,outfs,freqRange,getenvelope)
% d: recorded data, channel x sample
% fs: sampling rate of the data

% Nima, nimail@gmail.com

defaultfs = 400; % Hz

if ~exist('freqRange','var') || isempty(freqRange)
    freqRange=[70 150]; dispDefaultMessage(freqRange,'freqRange');
end
if ~exist('getenvelope','var') || isempty(getenvelope)
    getenvelope=1; dispDefaultMessage(getenvelope,'getenvelope');
end

if infs ~= defaultfs
    d = resample(d',defaultfs,infs)';
    fs = defaultfs;
end

% apply notch filter:
notchFreq=60;
while notchFreq<fs/2
    [b,a]=fir2(1000,[0 notchFreq-1 notchFreq-.5 notchFreq+.5 notchFreq+1 fs/2]/(fs/2),[1 1 0 0 1 1 ]);
    d=filtfilt(b,a,d')';
    notchFreq=notchFreq+60;
end
% calculate hilbert envelope:
if getenvelope
    [dh,cfs,sigma_fs] = CUprocessingHilbertTransform_filterbankGUI(d, fs, freqRange);
    dh2 = mean(abs(dh),3);
else
    [b,a]=butter(10,freqRange/(fs/2));
%     h = fvtool(b,a);
    dh2 = filtfilt(b,a,d);
end

dh2 = resample(dh2',outfs,fs)';
