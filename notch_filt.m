function signal = notch_filt(signal, P, sr)

% Uses 60 Hz noise to detect good/bad channels
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

n_channels = size(signal,2);

% notch filter parameters
b = cell(1,length(P.notch_freqs)); 
a = cell(1,length(P.notch_freqs));
for i = 1:length(P.notch_freqs)
    [b{i},a{i}] = iirnotch(P.notch_freqs(i) / (sr/2), P.notch_bw / (sr/2));
end

% fvtool(b{2}, a{2}, 'Fs', ecog_sr);
% figure;
% [h,t] = impz(b{1},a{1},[],ecog_sr);
% plot(t,h);
% xlim([0 1])

% apply notch filter
for i = 1:n_channels,
    for j = 1:length(P.notch_freqs)
        signal(:,i) = filtfilt(...
            b{j},a{j},signal(:,i));
    end
end
