close all; clear all; clc

% Nima's trick:  plot(10*tanh(b/10))

%% Params
subject = '053_NY625';
task = 'Statistical'; % only necessary if saving to HTK
block = 1;
dataorhtk = 'htk';
fs = 512;
keepallfigs = 0; % if 1, figures aren't automatically closed
stopandcheck = 1; % if 1, stops for user inspection after every cleaning type


%% Which parts to try

loaddata = 1;
assessaudio = 1;
assess


%% Set up cleanvars struct
cleanvars(1).bestaud = [];
cleanvars.filters = [];
cleanvars.channels = [];
cleanvars.amplifiers = [];
cleanvars.reftype = [];
cleanvars.epochs = [];
cleanvars.normdata = [];

%% Load data

switch dataorhtk
    case 'data'
        dpath = ['/Users/LauraLong/Documents/Lab/ECoG Data/' subject filesep 'original' filesep 'B' num2str(block)];
        load([dpath filesep 'data.mat'])
    case 'htk' 
        dpath = ['/Users/LauraLong/Documents/Lab/ECoG Data/' subject filesep 'processed' filesep 'B' num2str(block) ];
        channels = getchannelinds([dpath filesep 'htkraw'], 'htk');
        data = [];
        for k = 1:length(channels)
            [data(k,:),fs_data] = readhtk([dpath filesep 'htkraw' filesep 'Ch' num2str(channels(k)) '.htk']);
        end
        audchannels = getchannelinds([dpath filesep 'analog'], 'htk');
        aud = [];
        for k = 1:length(audchannels)
            [aud(k,:),fs_aud] = readhtk([dpath filesep 'analog' filesep 'a' num2str(audchannels(k)) '.htk']);
        end
        
end

%% View audio channels

if stopandcheck || keepallfigs
    figure;
    for i = 1:size(aud,1)
        subplot(size(aud,1),1,i)
        plot(aud(i,:)); title(num2str(i));
    end
    suptitle('Audio Channels');
end

%% Select best audio channel; view all data

bestaud = 11;


if stopandcheck || keepallfigs
    figure;
    subplot(211)
    imagesc(data); title('Raw Data');
    subplot(212);
    plot(aud(bestaud,:)); title('Audio');
    suptitle('Data and Best Audio');
    colormap jet
end

data_orig = data;

if stopandcheck, disp('Check audio results before continuing.'); keyboard, end
if ~keepallfigs, close all, end
aud = aud(bestaud,:);
cleanvars.bestaud = bestaud;

%% REMOVE DC 
%% Design high pass filter

highpassfreq = 0.5;

[b,a] = butter(1,highpassfreq/(fs/2),'high');
data_filtered = filtfilt(b,a,data')';

if stopandcheck || keepallfigs
    figure;
    subplot(211)
    plot(data_orig'); title('Unfiltered Data');
    subplot(212)
    plot(data_filtered'); title(['High Pass Filter at ' num2str(highpassfreq) ' Hz'])
    figure;
    freqz(b,a); title('High Pass Filter');
end

%% Look at filtered data

highpass = 0;

if highpass
    data = data_filtered; 
    cleanvars.filters.highpass = true;
    cleanvars.filters.highpassfreq = highpassfreq;
    cleanvars.filters.highpass_a = a;
    cleanvars.filters.highpass_b = b;
else
    cleanvars.filters.highpass = false;
    cleanvars.filters.highpassfreq = [];
end

if stopandcheck || keepallfigs
    figure;
    subplot(211)
    imagesc(data_orig); title('Original'); colorbar
    subplot(212)
    imagesc(data); title('After Filtering'); colorbar
end

if stopandcheck, disp('Check filtering results before continuing.'); keyboard, end
if ~keepallfigs, close all, end


%% CHANNELS
%% View single channels

if stopandcheck || keepallfigs
    channels = [101:110];
    for i = 1:length(channels)
        channel = channels(i);
        if highpass % in this case, show channels before and after filtering
            figure;
            subplot(221);
            plot(data_orig(channel,:)); title('Before Filtering: Waveform');
            subplot(223);
            pwelch(data_orig(channel,:));
            suptitle(['Channel ' num2str(channel)]);
            subplot(222);
            plot(data(channel,:)); title('After Filtering: Waveform');
            subplot(224);
            pwelch(data(channel,:));
        else
            figure;
            subplot(211);
            plot(data(channel,:)); title('Waveform');
            subplot(212);
            pwelch(data(channel,:));
            suptitle(['Channel ' num2str(channel)]);
        end
    end
end

%% Try removing channels

badchannels = [];
iffychannels = [46 47 54 70 95 105:108];
cleanestchannels = [];
data_cleanchannels = data;
data_cleanchannels(badchannels,:) = [];
data_bestchannels = data;
data_bestchannels(cat(2,badchannels,iffychannels),:) = [];

% keep the original indices
origchannelinds = 1:size(data,1);
origchannelinds(badchannels) = [];

if stopandcheck || keepallfigs
    figure;
    imagesc(data_cleanchannels); title('Clean Channels');
    colormap jet
    
    figure;
    subplot(131)
    imagesc(data_bestchannels); title('Good Channels');
    subplot(132);
    imagesc(data(badchannels,:)); title('Bad Channels');
    subplot(133);
    imagesc(data(iffychannels,:)); title('Iffy Channels');
    colormap jet
end

%% View all clean channel data

if stopandcheck || keepallfigs
    figure;
    subplot(311)
    plot(data_cleanchannels'); title('Waveform');
    subplot(312)
    imagesc(data_cleanchannels); title('Imagesc');
    subplot(313)
    plot(aud'); title('Audio');
    suptitle('Clean Channel Data');
    colormap jet
end

if stopandcheck, disp('Check channel results before continuing.'); keyboard, end
if ~keepallfigs, close all, end
cleanvars.channels.bad = badchannels;
cleanvars.channels.iffy = iffychannels;
cleanvars.channels.best = cleanestchannels;

%% REFERENCING
%% Try different referencing schemes

amplifiers = {1:32 33:64 65:96 97:128 129:160 161:192 193:224 225:256}; % where are the amplifiers split
amplifiers = {1:100 101:126};
amplifiers = {1:64 65:72 73:79 80:83 84:87 88:91 92:95 96:99 100:105 106:110 111:118 119:126};
amplifiers = {1:60 61:64 65:68 69:72 73:76 77:84 85:92 93:100 101:108 109:116 117:126};
amplifiers = {1:64 65:68 69:72 73:76 77:80 81:86 87:92 93:100 101:108 109:114 115:118 119:124 125:126};
amplifiers = {1:8 9:18 19:30 31:42 43:50 51:60 61:70 71:82};
amplifiers = {1:8 9:16 17:24 25:29 30:34 35:43 44:53 54:63 64:75 76:81}; 
amplifiers = {1:64 65:68 69:72 73:76 77:84 85:90 91:98 99:102 103:110 111:116 117:120 121:126 127:134 135:140 141:142 143:146 147:148};
amplifiers = {1:6 7:10 11:14 15:22 23:30 31:36 37:42 43:48 49:54 55:58 59:62 63:68 69:74 75:80 81:88 89:92 93:100 101:108 109:116 117:124 125:126};
amplifiers = {1:8 9:16 17:24 25:34 35:44 45:52 53:62 63:72 73:77 78:85 86:91};
amplifiers = {1:64 65:70 71:76 77:82 83:88 89:94};
amplifiers = {1:55 56:61 62:65 66:69 70:73 74:79 80:87 88:95 96:100 101:104 105:110};
maxPCs = 1; % how many PCs do you want to try

% Fixing amplifiers
origamplifiers = amplifiers;
for m = 1:length(badchannels)
    for n = 1:length(amplifiers)
        foundbadchannel = amplifiers{n} == badchannels(m);
        amplifiers{n}(foundbadchannel) = [];
        if isempty(amplifiers{n})
            amplifiers(n) = [];
        end
    end
end
amplifiersbyorigchannel = amplifiers;
currind = 0;
for m = 1:length(amplifiers)
    amplifiers{m} = [1:length(amplifiers{m})] + currind;
    try currind = amplifiers{m}(end); end
end

% ALL
% PCA (one component)
datastd = mapstd(data_cleanchannels);
[u,~] = svd(datastd*datastd');
dataref_pca_all = zeros(size(data_cleanchannels,1),size(data_cleanchannels,2),maxPCs);
for i = 1:maxPCs
    dataref_pca_all(:,:,i) = datastd - u(:,1:i)*(datastd'*u(:,1:i))';
end
% Average
dataref_avg_all = data_cleanchannels - repmat(mean(data_cleanchannels,1),[size(data_cleanchannels,1) 1]);
% Median
dataref_med_all = data_cleanchannels - repmat(median(data_cleanchannels,1),[size(data_cleanchannels,1) 1]);

if stopandcheck || keepallfigs
    figure;
    subplot(411);
    imagesc(data_cleanchannels); title('No Reference');
    subplot(412);
    imagesc(dataref_pca_all(:,:,1));  title(['PCA: First Component']);
    subplot(413);
    imagesc(dataref_avg_all); title('Average');
    subplot(414);
    imagesc(dataref_med_all); title('Median');
    suptitle('Reference Type: All');
    colormap jet
    if maxPCs > 1
        figure;
        subplot(maxPCs+1,1,1);
        imagesc(data_cleanchannels); title('original');
        for i = 1:maxPCs
            subplot(maxPCs+1,1,i+1)
            imagesc(dataref_pca_all(:,:,i)); title([num2str(i) ' components']);
        end
        suptitle('PCA Reference Comparing # Components: RefType All');
        colormap jet
    end
end

% BY AMPLIFIER
% PCA (one component)
datastd = mapstd(data_cleanchannels);
dataref_pca_amp = zeros(size(data_cleanchannels,1), size(data_cleanchannels,2),maxPCs);
for j = 1:length(amplifiers)
    datalocal = datastd(amplifiers{j},:);
    [u,~] = svd(datalocal*datalocal');
    for k = 1:maxPCs
        dataref_pca_amp(amplifiers{j},:,k) = datalocal - u(:,1:k) * (datalocal'*u(:,1:k))';
    end
end
% Average
dataref_avg_amp = zeros(size(data_cleanchannels));
for j = 1:length(amplifiers)
    datalocal = data_cleanchannels(amplifiers{j},:);
    dataref_avg_amp(amplifiers{j},:) = datalocal - repmat(mean(datalocal,1),[size(datalocal,1) 1]);
end
% Median
dataref_med_amp = zeros(size(data_cleanchannels));
for j = 1:length(amplifiers)
    datalocal = data_cleanchannels(amplifiers{j},:);
    dataref_med_amp(amplifiers{j},:) = datalocal - repmat(median(datalocal,1),[size(datalocal,1) 1]);
end

if stopandcheck || keepallfigs
    figure;
    subplot(411);
    imagesc(data_cleanchannels); title('No Reference');
    subplot(412);
    imagesc(dataref_pca_amp(:,:,1));  title(['PCA: First Component']);
    subplot(413);
    imagesc(dataref_avg_amp); title('Average');
    subplot(414);
    imagesc(dataref_med_amp); title('Median');
    suptitle('Reference Type: Amp');
    colormap jet
    if maxPCs > 1
        figure;
        subplot(maxPCs+1,1,1);
        imagesc(data_cleanchannels); title('original');
        for i = 1:maxPCs
            subplot(maxPCs+1,1,i+1)
            imagesc(dataref_pca_amp(:,:,i)); title([num2str(i) ' components']);
        end
        suptitle('PCA Reference Comparing # Components: RefType Amp');
        colormap jet
    end
end

% LOCALLY
numchannels = size(data_cleanchannels,1);
% Average
dataref_avg_loc = zeros(size(data_cleanchannels));
for j = 1:numchannels
    thisdata = data_cleanchannels(j,:);
    switch j
        case 1
            datalocal = data_cleanchannels(j+1,:);
        case numchannels
            datalocal = data_cleanchannels(j-1,:);
        otherwise
            datalocal = data_cleanchannels([j-1 j+1],:);
    end
    dataref_avg_loc(j,:) = thisdata - mean(datalocal,1);
end
% Median
dataref_med_loc = zeros(size(data_cleanchannels));
for j = 1:numchannels
    thisdata = data_cleanchannels(j,:);
    switch j
        case 1
            datalocal = data_cleanchannels(j+1,:);
        case numchannels
            datalocal = data_cleanchannels(j-1,:);
        otherwise
            datalocal = data_cleanchannels([j-1 j+1],:);
    end
    dataref_med_loc(j,:) = thisdata - mean(datalocal,1);
end

if stopandcheck || keepallfigs
    figure;
    subplot(311);
    imagesc(data_cleanchannels); title('No Reference');
    subplot(312);
    imagesc(dataref_avg_loc); title('Average');
    subplot(313);
    imagesc(dataref_med_loc); title('Median');
    suptitle('Reference Type: Local');
    colormap jet
end

%% Select best referencing

bestreftype = 'amp'; % can be all, amp, local, or none
bestrefsubtype = 'med'; % avg, pca, or med
bestPC = 2;

switch lower(bestreftype(1:3))
    case 'all'
        switch lower(bestrefsubtype)
            case 'pca'
                dataref_pca_all = dataref_pca_all(:,:,bestPC);
                data_referenced = dataref_pca_all;
            case 'avg'
                data_referenced = dataref_avg_all;
            case 'med'
                data_referenced = dataref_med_all;
        end
    case 'amp'
        switch lower(bestrefsubtype)
            case 'pca'
                dataref_pca_amp = dataref_pca_amp(:,:,bestPC);
                data_referenced = dataref_pca_amp;
            case 'avg'
                data_referenced = dataref_avg_amp;
            case 'med'
                data_referenced = dataref_med_amp;
        end
    case 'loc'
        switch lower(bestrefsubtype)
            case 'avg'
                data_referenced = dataref_avg_loc;
            case 'med'
                data_referenced = dataref_med_loc;
        end
    case 'non'
        data_referenced = data_cleanchannels;
end

if stopandcheck || keepallfigs
    figure;
    subplot(211);
    imagesc(data_cleanchannels); title('Unreferenced');
    subplot(212);
    imagesc(data_referenced); title(['Reference Type: ' bestreftype ' ' bestrefsubtype]);
    colormap jet
end

if stopandcheck, disp('Check referencing results before continuing.'); keyboard, end
if ~keepallfigs, close all, end
cleanvars.reftype = [bestreftype(1:3) ' ' bestrefsubtype];
if strcmpi(bestrefsubtype,'pca')
    cleanvars.reftype = [cleanvars.reftype num2str(bestPC)];
end
cleanvars.amplifiers = origamplifiers;

%% EPOCHS
%% Calculate mean and PC for various epochs

epochsize = 1000;

refepochs = [1:epochsize:size(data_referenced,2) size(data_referenced,2)];
epochmeans = zeros(size(data_referenced,1),length(refepochs)-1);
epochstds = zeros(size(epochmeans));
epochmeans_preref = zeros(size(epochmeans));
epochstds_preref = zeros(size(epochmeans));
for i = 1:length(refepochs)-1
    thisepoch = refepochs(i):refepochs(i+1);
    epochmeans(:,i) = mean(data_referenced(:,thisepoch),2);
    epochstds(:,i) = std(data_referenced(:,thisepoch),[],2);
    
    epochmeans_preref(:,i) = mean(data_cleanchannels(:,thisepoch),2);
    epochstds_preref(:,i) = std(data_cleanchannels(:,thisepoch),[],2);
end

if stopandcheck || keepallfigs
    figure;
    subplot(221)
    imagesc(epochmeans); title('Mean by Epoch: PostRef');
    subplot(223)
    imagesc(epochstds); title('STD by Epoch: PostRef');
    subplot(222)
    imagesc(epochmeans_preref); title('Mean by Epoch: PreRef');
    subplot(224)
    imagesc(epochstds_preref); title('STD by Epoch: PreRef');
    colormap jet
end


%% View single epochs
if stopandcheck || keepallfigs
    epochs = {1:10000};
    for i = 1:length(epochs)
        epoch = epochs{i};
        figure;
        subplot(311);
        plot(data_referenced(:,epoch)'); title('Waveform');
        subplot(312);
        imagesc(data_referenced(:,epoch)); title('Imagesc');
        subplot(313);
        plot(aud(:,epoch)'); title('Audio');
        suptitle(['Epoch: samples ' num2str(epoch(1)) ' to ' num2str(epoch(end))]);
        colormap jet
    end
end

%% Try removing large artifacts

badepochs = [];
data_cleanepochs = data_referenced;
data_cleanepochs(:,badepochs) = [];
aud_cleanepochs = aud;
aud_cleanepochs(:,badepochs) = [];

if stopandcheck || keepallfigs
    figure;
    subplot(311)
    plot(data_cleanepochs'); xlim([0 size(data_cleanepochs,2)]); title('Waveform');
    subplot(312)
    imagesc(data_cleanepochs); title('Imagesc');
    subplot(313)
    plot(aud_cleanepochs'); xlim([0 size(aud_cleanepochs,2)]); title('Audio');
    suptitle('Clean Epochs');
    colormap jet
end

if stopandcheck, disp('Check epoch results before continuing.'); keyboard, end
if ~keepallfigs, close all, end

%% NORMALIZE
%% View z-scored and non z-scored data

data_norm = zscore(data_cleanepochs,[],2);
if stopandcheck || keepallfigs
    figure;
    subplot(221)
    plot(data_cleanepochs'); title('Plain');
    subplot(223)
    imagesc(data_cleanepochs);
    subplot(222)
    plot(data_norm'); title('Z-Scored');
    subplot(224)
    imagesc(data_norm);
    suptitle('Normalized Data');
    colormap jet
end
%% Decide whether to z-score

normdata = 0;

if normdata
    cleanvars.normdata = true;
    data_final = data_norm;
    if stopandcheck || keepallfigs
    figure;
    subplot(311)
    plot(data_cleanchannels'); title('Waveform');
    subplot(312)
    imagesc(data_cleanchannels); title('Imagesc');
    subplot(313)
    plot(aud'); title('Audio');
    suptitle('Normalized Data');
    colormap jet
    end
else
    cleanvars.normdata = false;
    data_final = data_cleanepochs;
end

if stopandcheck, disp('Check normalization results before continuing.'); keyboard, end
if ~keepallfigs, close all, end
cleanvars.normdata = normdata;


%% Display each cleaning step

figure;
subplot(611);
imagesc(data_orig); title('original');
subplot(612)
imagesc(data); title('filtered');
subplot(613)
imagesc(data_cleanchannels); title('without bad channels');
subplot(614)
imagesc(data_referenced); title(['after ' bestreftype ' ' bestrefsubtype ' referencing']);
subplot(615)
imagesc(data_cleanepochs); title('without bad epochs');
subplot(616)
imagesc(data_final); title('after normalization');
suptitle('Results of Cleaning Data');
colormap jet

if stopandcheck, keyboard, disp('Check overall results before continuing.'); end

%% Save results and make HTK if desired

saveinput = input('Save found params? ');
if saveinput == 1

    save(['/Users/LauraLong/Documents/Lab/ECoG Data/' subject filesep 'processed' filesep 'B' num2str(block) filesep subject '_dataqualityvars.mat'],'cleanvars','dataorhtk','subject','block');
    
    savedata = input('Save cleaned data? ');
    if savedata == 1
        save(['/Users/LauraLong/Documents/Lab/ECoG Data/' subject filesep 'processed' filesep 'B' num2str(block) filesep subject '_cleaneddata.mat'],'data_cleanepochs','subject','block');
    end
    
    htkinput = input('Make HTK based on these params? ');
    if htkinput == 1
        dpath = ['/Users/LauraLong/Documents/Lab/ECoG Data/' subject];
        soundpath = ['/Users/LauraLong/Documents/Lab/ECoG Data/Task_Sounds/' task];
        channels = [];
        whichdata = [];
        dataf = [];
        HTKtoCleanData(dpath,whichdata,block,channels,dataf,cleanvars.filters.highpassfreq,cleanvars.reftype,cleanvars.amplifiers,cleanvars.channels.bad,cleanvars.epochs,cleanvars.channels.iffy,cleanvars.normdata,1);
    end
    
end
