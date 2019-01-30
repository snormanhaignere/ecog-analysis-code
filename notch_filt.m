function signal = notch_filt(signal, sr, varargin)

% Uses 60 Hz noise to detect good/bad channels
% 
% Derived from Schalk lab code
% 
% signal: [time x electrode] singal matrix
% 
% sr: sampling rate of the signal
% 
% 2016-08-12: Created, Sam NH
% 
% 2019-01-21: Altered parameter handling, Sam NH

I.bw = 1;
I.freqs = [60, 120, 160, 240];
I = parse_optInputs_keyvalue(varargin, I);

n_channels = size(signal,2);

% notch filter parameters
b = cell(1,length(I.freqs)); 
a = cell(1,length(I.freqs));
for i = 1:length(I.freqs)
    [b{i},a{i}] = iirnotch(I.freqs(i) / (sr/2), I.bw / (sr/2));
end

% fvtool(b{2}, a{2}, 'Fs', ecog_sr);
% figure;
% [h,t] = impz(b{1},a{1},[],ecog_sr);
% plot(t,h);
% xlim([0 1])

% apply notch filter
for i = 1:n_channels
    for j = 1:length(I.freqs)
        signal(:,i) = filtfilt(...
            b{j},a{j},signal(:,i));
    end
end
