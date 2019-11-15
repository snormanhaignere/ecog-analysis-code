function [ cleanvars ] = findCleanDataVariables(datapath,block,task)
% Runs through what ECoGdataqualitytools used to do

fprintf('\nfindCleanDataVariables\n');


%% Load data

datamatpath = [datapath filesep 'original' filesep 'B' num2str(block) filesep 'data.mat'];
if exist(datamatpath,'file')
    dataorhtk = input('Enter 1 to load from data.mat and 2 to load from htk: ');
else
    dataorhtk = 2;
end

switch dataorhtk
    case 1
        load(datamatpath);
        data_orig = data;
    case 2 
        dpath = [datapath filesep 'processed' filesep 'B' num2str(block) ];
        channels = getchannelinds([dpath filesep 'htkraw'], 'htk');
        data_orig = [];
        for k = 1:length(channels)
            [data_orig(k,:),fs_data] = readhtk([dpath filesep 'htkraw' filesep 'Ch' num2str(channels(k)) '.htk']);
        end
        data = data_orig;
        audchannels = getchannelinds([dpath filesep 'analog'], 'htk');
        aud = [];
        for k = 1:length(audchannels)
            [aud(k,:),fs_aud] = readhtk([dpath filesep 'analog' filesep 'a' num2str(audchannels(k)) '.htk']);
        end
        
end

numversions = 1;

%% Audio 

runsection = input('Visualize audio to select best channel? '); 

if runsection
    % Plot every audio channel
    figure;
    for i = 1:size(aud,1)
        subplot(size(aud,1),1,i)
        plot(aud(i,:)); title(num2str(i));
    end
    suptitle('Audio Channels');
    
    % Find best audio
    bestaud = input('  Enter best audio channel: ');
    
    % Now replot with best audio
    figure;
    subplot(211)
    imagesc(data); title('Raw Data');
    subplot(212);
    plot(aud(bestaud,:)); title('Audio');
    suptitle('Data and Best Audio');
    colormap jet
    
else % otherwise just select the first channel
    bestaud = 1;
end

% Keep best audio and update cleanvars
aud = aud(bestaud,:);
cleanvars.bestaud = bestaud;

%% Remove DC

runsection = input('Try high pass filter to remove DC drift? '); close all; 

if runsection
    
    % Design 0.5Hz high pass filter
    highpassfreq = input('   Enter high pass frequency: ');
%     highpassfreq = 0.5;
    [b,a] = butter(1,highpassfreq/(fs_data/2),'high');
    data_filtered = filtfilt(b,a,data')';
    
    % Plot filtered and unfiltered data
    figure;
    subplot(211)
    plot(data'); title('Unfiltered Data');
    subplot(212)
    plot(data_filtered'); title(['High Pass Filter at ' num2str(highpassfreq) ' Hz'])
    figure;
    freqz(b,a); title('High Pass Filter');
    
    % Decide whether to keep audio
    highpass = input('  Filter data? ');
    
    % If filtering, save choices and visualize
    if highpass
        data = data_filtered;
        cleanvars.filters.highpass = true;
        cleanvars.filters.highpassfreq = highpassfreq;
        cleanvars.filters.highpass_a = a;
        cleanvars.filters.highpass_b = b;
        
        %Plot
        figure;
        subplot(211)
        imagesc(data_orig); title('Original'); colorbar
        subplot(212)
        imagesc(data); title('After Filtering'); colorbar
        
        numversions = numversions+1;
        
    else % otherwise keep current data
        data = data_orig;
        cleanvars.filters.highpass = false;
        cleanvars.filters.highpassfreq = [];
    end
    
else
    data = data_orig;
    cleanvars.filters.highpass = false;
    cleanvars.filters.highpassfreq = [];
end


%% Channel selection

runsection = input('Look for noisy channels? '); close all; 

if runsection
    
    % Disp number of channels, set checkchannels to 1 and bad/iffy channels to empty
    disp(['  Total number of data channels: ' num2str(size(data,1))]);
    checkchannels = 1;
    badchannels = [];
    iffychannels = [];
    
    % While the user wants to continue, plot chunks of channels and ask which are iffy/bad and add to lists
    while checkchannels == 1
        if ~exist('chanstoplot','var') || isequal(chanstoplot,1);
            chanstoplot = input('  Enter channels to plot: ');
        end
        for i = 1:length(chanstoplot)
            chantoplot = chanstoplot(i);
            if cleanvars.filters.highpass % in this case, show channels before and after filtering
                figure;
                subplot(221);
                plot(data_orig(chantoplot,:)); title('Before Filtering: Waveform');
                subplot(223);
                pwelch(data_orig(chantoplot,:));
                suptitle(['Channel ' num2str(chantoplot)]);
                subplot(222);
                plot(data(chantoplot,:)); title('After Filtering: Waveform');
                subplot(224);
                pwelch(data(chantoplot,:));
            else
                figure;
                subplot(211);
                plot(data(chantoplot,:)); title('Waveform');
                subplot(212);
                pwelch(data(chantoplot,:));
                suptitle(['Channel ' num2str(chantoplot)]);
            end
        end
        
        newbadchan = input('   Enter numbers of bad channels: ');
        badchannels = [badchannels newbadchan];
        newiffychan = input('   Enter numbers of iffy channels: ');
        iffychannels = [iffychannels newiffychan];
        
        checkchannels = input('  Check additional channels? ');
        if ~isempty(checkchannels) && ~isequal(checkchannels,0);
            chanstoplot = checkchannels;
            checkchannels = 1;
        else
            checkchannels = 0;
        end 
        
    end
    
    % Once complete, delete bad channels from data
    data_cleanchannels = data;
    data_cleanchannels(badchannels,:) = [];
    origchannelinds = 1:size(data,1);
    origchannelinds(badchannels) = [];
    
    % Plot clean data
    figure;
    imagesc(data_cleanchannels); title('Clean Channels');
    colormap jet
    figure;
    subplot(131)
    imagesc(data_cleanchannels); title('Good Channels');
    subplot(132);
    imagesc(data(badchannels,:)); title('Bad Channels');
    subplot(133);
    imagesc(data(iffychannels,:)); title('Iffy Channels');
    colormap jet
    
    % Save data
    data = data_cleanchannels;
    numversions = numversions+1;
    
else % otherwise
    origchannelinds = 1:size(data,1);
    badchannels = [];
    iffychannels = [];
    
end

badchannels = sort(unique(badchannels),'ascend');
iffychannels = sort(unique(iffychannels),'ascend');


% Save information into cleanvars
cleanvars.channels.bad = badchannels;
cleanvars.channels.iffy = iffychannels;


%% Referencing Schemes 

runsection = input('Try different referencing schemes? '); close all;

if runsection
    
    % Determine whether to do PCs
    maxPCs = input('  Enter max # of PCs to try: ');
    
    % Global Referencing
    runsubsection = input('  Try global referencing? ');
    if runsubsection
        
        % Find results
        if maxPCs>0
            % PCA (one component)
            datastd = mapstd(data);
            [u,~] = svd(datastd*datastd');
            dataref_pca_all = zeros(size(data,1),size(data,2),maxPCs);
            for i = 1:maxPCs
                dataref_pca_all(:,:,i) = datastd - u(:,1:i)*(datastd'*u(:,1:i))';
            end
        end
        % Average
        dataref_avg_all = data - repmat(mean(data,1),[size(data,1) 1]);
        % Median
        dataref_med_all = data - repmat(median(data,1),[size(data,1) 1]);
        
        % Plot
        figure;
        subplot(411);
        imagesc(data); title('No Reference');
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
            imagesc(data); title('original');
            for i = 1:maxPCs
                subplot(maxPCs+1,1,i+1)
                imagesc(dataref_pca_all(:,:,i)); title([num2str(i) ' components']);
            end
            suptitle('PCA Reference Comparing # Components: RefType All');
            colormap jet
        end
    end
    
    % Amplifier Referencing
    runsubsection = input('  Try amplifier referencing? ');
    if runsubsection
    
        % Need to get amplifier information- first load and display chnames
        ampsknown = input('   Do you know the amplifier indices? ');
        if isequal(ampsknown,0)
            chnamefile = [datapath filesep 'original' filesep 'B' num2str(block) filesep 'allchnames_B' num2str(block) '.mat'];
            load(chnamefile)
            disp(chnames');
            disp('   Determine amplifiers from chnames above before continuing.'); keyboard;
            amplifiers = input('   Enter amplifier indices: ');
        else
            amplifiers = ampsknown;
        end
        
        
        % Now determine proper amplifier information based on badchannels
        % Delete bad channels from each amplifier
        origamplifiers = amplifiers;
        for m = 1:length(badchannels)
            for n = 1:length(origamplifiers)
                foundbadchannel = amplifiers{n} == badchannels(m);
                amplifiers{n}(foundbadchannel) = [];
            end
        end
        % Delete empty amplifiers from amplifier list
        foundemptyamp = false(1,length(amplifiers));
        for m = 1:length(amplifiers)
            foundemptyamp(m) = isempty(amplifiers{m});
        end
        amplifiers(foundemptyamp) = [];
        % Renumber amplifier channels from 1:end rather than having gaps
        amplifiersbyorigchannel = amplifiers;
        currind = 0;
        for m = 1:length(amplifiers)
            amplifiers{m} = [1:length(amplifiers{m})] + currind;
            try currind = amplifiers{m}(end); end
        end
        
        % Find results
        % PCA
        if maxPCs>0
            datastd = mapstd(data);
            dataref_pca_amp = zeros(size(data,1), size(data,2),maxPCs);
            for j = 1:length(amplifiers)
                datalocal = datastd(amplifiers{j},:);
                [u,~] = svd(datalocal*datalocal');
                for k = 1:maxPCs
                    dataref_pca_amp(amplifiers{j},:,k) = datalocal - u(:,1:k) * (datalocal'*u(:,1:k))';
                end
            end
        end
        % Average
        dataref_avg_amp = zeros(size(data));
        for j = 1:length(amplifiers)
            datalocal = data(amplifiers{j},:);
            dataref_avg_amp(amplifiers{j},:) = datalocal - repmat(mean(datalocal,1),[size(datalocal,1) 1]);
        end
        % Median
        dataref_med_amp = zeros(size(data));
        for j = 1:length(amplifiers)
            datalocal = data(amplifiers{j},:);
            dataref_med_amp(amplifiers{j},:) = datalocal - repmat(median(datalocal,1),[size(datalocal,1) 1]);
        end
        
        % Plot
        figure;
        subplot(411);
        imagesc(data); title('No Reference');
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
    
    
    % Ask which is best
    bestreftype = input('  Enter best reference type (all, amp, or none): ','s');
    bestrefsubtype = input('  Enter best reference subtype (avg, pca, or med): ','s');
    
    % Keep desired data
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
    
    % Plot 
    figure;
    subplot(211);
    imagesc(data_cleanchannels); title('Unreferenced');
    subplot(212);
    imagesc(data_referenced); title(['Reference Type: ' bestreftype ' ' bestrefsubtype]);
    colormap jet
    
    % Save data
    data = data_referenced;
    numversions = numversions+1;

else
    bestreftype = 'none';
    bestrefsubtype = [];
    origamplifiers = [];
end

cleanvars.reftype = [bestreftype(1:3) ' ' bestrefsubtype];
if strcmpi(bestrefsubtype,'pca')
    cleanvars.reftype = [cleanvars.reftype num2str(bestPC)];
end
if exist('origamplifiers','var')
    cleanvars.amplifiers = origamplifiers;
else
    cleanvars.amplifiers = [];
end


%% Epochs 

runsection = input('Check epochs for noise? ');

if runsection
    
    % First show statistics by epoch; allow user to determine size
    epochsize = input('  Enter epoch size to check stats: ');
    refepochs = [1:epochsize:size(data,2) size(data,2)];
    epochmeans = zeros(size(data,1),length(refepochs)-1);
    epochstds = zeros(size(epochmeans));
    for i = 1:length(refepochs)-1
        thisepoch = refepochs(i):refepochs(i+1);
        epochmeans(:,i) = mean(data(:,thisepoch),2);
        epochstds(:,i) = std(data(:,thisepoch),[],2);
    end
    % Plot these results
    figure;
    subplot(211)
    imagesc(epochmeans); title('Mean by Epoch: PostRef');
    subplot(212)
    imagesc(epochstds); title('STD by Epoch: PostRef');
    colormap jet
    
    % Look at specific epochs
    checkepochs = input('  Check specific epochs? ');
    badepochs = [];
    while checkepochs
        thisepoch = input('   Enter epoch to check: ');
        figure;
        subplot(311);
        plot(data(:,thisepoch)'); title('Waveform');
        subplot(312);
        imagesc(data(:,thisepoch)); title('Imagesc');
        subplot(313);
        plot(aud(:,thisepoch)'); title('Audio');
        suptitle(['Epoch: samples ' num2str(thisepoch(1)) ' to ' num2str(thisepoch(end))]);
        colormap jet
        
        newbadepoch = input('   Enter samples to cut: ');
        badepochs = [badepochs newbadepoch];
        
        checkepochs = input('  Check additional epochs? ');
    end
    
    % Eliminate chosen epochs
    data_cleanepochs = data;
    data_cleanepochs(:,badepochs) = [];
    aud_cleanepochs = aud;
    aud_cleanepochs(:,badepochs) = [];
    
    % Plot
    figure;
    subplot(311)
    plot(data_cleanepochs'); xlim([0 size(data_cleanepochs,2)]); title('Waveform');
    subplot(312)
    imagesc(data_cleanepochs); title('Imagesc');
    subplot(313)
    plot(aud_cleanepochs'); xlim([0 size(aud_cleanepochs,2)]); title('Audio');
    suptitle('Clean Epochs');
    colormap jet
    
    
    % Save data
    data = data_cleanepochs;
    aud = aud_cleanepochs;
    numversions = numversions+1;
    
else
    badepochs = [];
    
end

cleanvars.epochs = badepochs;


%% Normalization

runsection = input('Try z-scoring by channel? '); close all;

if runsection

    % Try normalization
    data_norm = zscore(data,[],2);
    
    % Plot
    figure;
    subplot(221)
    plot(data'); title('Plain');
    subplot(223)
    imagesc(data);
    subplot(222)
    plot(data_norm'); title('Z-Scored');
    subplot(224)
    imagesc(data_norm);
    suptitle('Normalized Data');
    colormap jet
    
    % Ask whether to keep; if keep, save
    normdata = input('  Keep z-score? ');
    if normdata
        data = data_norm;
        cleanvars.normdata = true;
        numversions = numversions+1;
    else
        cleanvars.normdata = false;
    end
else
    cleanvars.normdata = false;
end


%% Plot all versions
close all; 

figure; 
subplot(numversions,1,1); imagesc(data_orig); title('original'); plotnum = 2;
if exist('data_filtered','var')
    subplot(numversions,1,plotnum); imagesc(data_filtered); title('filtered');
    plotnum = plotnum+1;
end
if exist('data_cleanchannels','var')
    subplot(numversions,1,plotnum); imagesc(data_cleanchannels); title('without bad channels');
    plotnum = plotnum+1;
end
if exist('data_referenced','var')
    subplot(numversions,1,plotnum); imagesc(data_referenced); title('after referencing');
    plotnum = plotnum+1;
end
if exist('data_cleanepochs','var')
    subplot(numversions,1,plotnum); imagesc(data_cleanepochs); title('without bad epochs');
    plotnum = plotnum+1;
end
if exist('data_norm','var')
    subplot(numversions,1,plotnum); imagesc(data_norm); title('after normalization');
    plotnum = plotnum+1;
end
suptitle('Results of Cleaning Data'); colormap jet



%% Save results if desired 
runsection = input('Save found params? ');
if runsection
    save([datapath filesep 'processed' filesep 'B' num2str(block) filesep 'dataqualityvars.mat'],'cleanvars','block');
    savedata = input('Save cleaned data? ');
    if savedata
        save(['/Users/LauraLong/Documents/Lab/ECoG Data/' subject filesep 'processed' filesep 'B' num2str(block) filesep subject '_cleaneddata.mat'],'data_cleanepochs','block');
    end
end

fprintf('Completed findCleanDataVariables\n');

end

