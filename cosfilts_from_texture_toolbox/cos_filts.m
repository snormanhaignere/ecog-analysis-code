function [filts,Hz_cutoffs,freqs] = ...
    cos_filts(signal_length, sr, Q, low_lim, hi_lim)

if rem(signal_length,2)==0 %even length
    nfreqs = signal_length/2;%does not include DC
    max_freq = sr/2;
    freqs = [0:max_freq/nfreqs:max_freq]; %go all the way to nyquist
else %odd length
    nfreqs = (signal_length-1)/2;
    max_freq = sr*(signal_length-1)/2/signal_length; %max freq is just under nyquist
    freqs = [0:max_freq/nfreqs:max_freq];
end   

log2(low_lim) : 

cos_filts = zeros(nfreqs+1,N);

if hi_lim>sr/2
    hi_lim = max_freq;
end

num_filters = 2*N+1;
%make cutoffs evenly spaced on an erb scale
spacing = (freq2erb(hi_lim)-freq2erb(low_lim))/(num_filters+1);%in ERBs
center_freqs = linspace(freq2erb(low_lim)+spacing, freq2erb(hi_lim)-spacing, num_filters); %in ERBs

for k=1:num_filters
    l = erb2freq(center_freqs(k)-2*spacing);
    h = erb2freq(center_freqs(k)+2*spacing);
    l_ind = min(find(freqs>l));
    h_ind = max(find(freqs<h));
    avg = (freq2erb(l)+freq2erb(h))/2;
    rnge = (freq2erb(h)-freq2erb(l));
    cos_filts(l_ind:h_ind,k) = cos((freq2erb( freqs(l_ind:h_ind) ) - avg)/rnge*pi); %map cutoffs to -pi/2, pi/2 interval
end

%add lowpass and highpass to get perfect reconstruction
filts = zeros(nfreqs+1,num_filters+4);
filts(:,3:num_filters+2) = cos_filts;
%lowpass filters go up to peaks of first, second cos filters
h_ind = max(find(freqs<erb2freq(center_freqs(1))));
filts(1:h_ind,1) = sqrt(1 - filts(1:h_ind,3).^2);
h_ind = max(find(freqs<erb2freq(center_freqs(2))));
filts(1:h_ind,2) = sqrt(1 - filts(1:h_ind,4).^2);
%highpass filters go down to peaks of last two cos filters
l_ind = min(find(freqs>erb2freq(center_freqs(num_filters))));
filts(l_ind:nfreqs+1,num_filters+4) = sqrt(1 - filts(l_ind:nfreqs+1,num_filters+2).^2);
l_ind = min(find(freqs>erb2freq(center_freqs(num_filters-1))));
filts(l_ind:nfreqs+1,num_filters+3) = sqrt(1 - filts(l_ind:nfreqs+1,num_filters+1).^2);

filts = filts/sqrt(2); %so that squared freq response adds to 1

center_freqs = erb2freq([center_freqs(1)-2*spacing center_freqs(2)-2*spacing center_freqs center_freqs(num_filters-1)+2*spacing center_freqs(num_filters)+2*spacing]);
Hz_cutoffs = center_freqs;
Hz_cutoffs(find(Hz_cutoffs<0)) = 1;

%subplot(2,1,1); plot(freqs,sum(filts.^2,2))
%subplot(2,1,2); semilogx(freqs,sum(filts.^2,2))
