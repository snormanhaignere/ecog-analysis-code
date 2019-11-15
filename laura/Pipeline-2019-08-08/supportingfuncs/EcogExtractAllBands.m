function dh2 = EcogExtractAllBands (d, infs, outfs)
% d: recorded data, channel x sample
% fs: sampling rate of the data

% Nima, nimail@gmail.com

defaultfs = 400; % Hz
freqRange=[1 150];

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
[dh,cfs,sigma_fs] = CUprocessingHilbertTransform_filterbankGUI(d, fs, freqRange);
%
bands = {1:5,6:11,12:19,20:26,27:34,35:42};
dh2 = [];
for cnt1 = 1:length(bands)
    dh2(:,cnt1) = mean(abs(dh(:,:,bands{cnt1})),3);
end
% dh2 = mean(abs(dh(:,:,[),3);
%dh2 = mapstd(dh2);
dh2 = resample(dh2,outfs,fs)';
