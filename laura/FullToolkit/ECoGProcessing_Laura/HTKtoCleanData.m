function CMRtag = HTKtoCleanData(basepath,whichdata,blocks,allchannels,dataf,highpassfreq,fullreftype,amplifiers,badchannels,badepochs,iffychannels,normdata,visualizeresults)
% function HTKtoCMR
%       Extracts high gamma from HTK files, writes to new HTK files
%
% Inputs:
%       basepath (optional): subject folder
%           basepath must contain htk data in a 'processed' subdirectory
%               eg: basepath / processed <- contains htkraw files
%           default: current directory
%       whichdata (optional): which processed file you want to use
%           default: htkraw
%       blocks (optional): specify which blocks to process
%           default: all blocks in basepath / processed
%       channels (optional): specify which channels to process
%           default: to all channels in each block
%       dataf (optional): specify data fs of source htk files
%           default: 1000Hz
%
% Result:
%       Creates 'highgamma' directory in each processed block file.
%       Extracts high gamma from files in basepath / processed / whichdata.
%       Writes a new .htk file with high gamma.
%
%       Resulting file structure:
%           From htkraw: basepath / processed / B# / highgamma / Ch#.htk
%           Data: basepath / processed / B# / whichdata_highgamma / Ch#.htk
%
% Hard-Coded Assumptions:
%       1. Default data channels are all channels in the first block
%       2. Original data is in htkraw folder.
%
%
% Written: 10/2015 by Laura
%
% Last Updated: 07/29/16 by Laura to add median option
%
%% Set inputs, find path

fprintf('\nHTKtoCleanData\n');

% Check if basepath was entered; if not, default to current directory
if ~exist('basepath','var') || isempty(basepath)
    basepath = pwd; % if not specified, default to the current directory
end

% Define datapath from basepath
datapath = [basepath filesep 'processed' filesep];

% Throw error if processed file isn't already there
if ~exist(datapath,'dir')
    error('Cannot find processed directory; check basepath');
end

% Check if whichdata was entered; if not, default to htkraw
if ~exist('whichdata','var') || isempty(whichdata)
    whichdata = 'htkraw'; % if not specified, default to the current directory
end

% Define appropriate directory tag
if strcmp(whichdata,'htkraw')
    dir_tag = []; % no tag if converting from raw
else
    dir_tag = [whichdata '_']; % include tag if using already processed data
end

% Set parameter defaults

% Data sampling rate defaults to 1000 (aud is hard coded to 2400)
if ~exist('dataf','var') || isempty(dataf) % sampling frequency default to 1000
    dataf = 1000;
end

% Blocks defaults to all blocks found in the original folder
if ~exist('blocks','var') || isempty(blocks) % blocks default to all
    blocks = getblockinds(datapath); %getblockinds finds vector of block numbers from basepath/processed/
end

% Channels defaults to all channels found in the first block folder
if ~exist('allchannels','var') || isempty(allchannels) % channels default to all
    allchannels = getchannelinds([datapath 'B' num2str(blocks(1)) filesep 'htkraw'],'htk');
end


if ~exist('badchannels','var') || isempty(badchannels)
    badchannels = [];
end
if ~exist('iffychannels','var') || isempty(iffychannels)
    iffychannels = [];
end

if ~exist('normdata','var') || isempty(normdata)
    normdata = 0;
end

if ~exist('visualizeresults','var') || isempty(visualizeresults)
    visualizeresults = 0;
end

if ~exist('amplifiers','var') || isempty(amplifiers)
    amplifiers = {allchannels};
end

if ~exist('highpassfreq','var') || isempty(highpassfreq)
    highpass = 0;
    highpassfreq = [];
else
    highpass = 1;
end

warning('Currently there''s a channel bug if run for more than one block at a time!');


% Set up reference types
if ~exist('fullreftype','var') || isempty(fullreftype) % if there's nothing there, default to none
    reftype = '';
    refsubtype = '';
else % otherwise parse into the main reftype and the subtype
    reftypes = strsplit(fullreftype,' '); % split by space bar
    reftype = lower(reftypes{1}); % main (all, local, or amps)
    refsubtype = lower(reftypes{2}); % subtype (pca, avg, med)
    if length(refsubtype) == 4 && strcmp(refsubtype(1:3),'pca')
        try
            bestPC = int8(str2num(refsubtype(4)));
        catch
            bestPC = 1;
        end
        refsubtype = 'pca';
    elseif strcmp(refsubtype(1:3),'pca')
        bestPC = 1;
    end
end


% Set up saved variable
cleanvars(1).filters = [];
cleanvars.filters.highpass = highpass;
cleanvars.channels.bad = badchannels;
cleanvars.channels.iffy = iffychannels;
cleanvars.amplifiers.allchannels = amplifiers;
cleanvars.reftype = fullreftype;
cleanvars.epochs = badepochs;
cleanvars.normdata = normdata;


%% Loop over blocks, creating CMR HTK files

for i = blocks
    
    % Generate block name
    blockname = ['B' int2str(i)]; % create block string for use in files
    disp(['Processing ' blockname]);
    
    % Throw error if whichdata file isn't already there
    if ~exist([datapath blockname, filesep, whichdata],'dir')
        error(['Cannot find ' whichdata ' directory; check processed folder']);
    end
    
    %% Load data
    
    % Load data
    disp('...loading data');
    data = [];
    for j = 1:length(allchannels)
        [data(j,:),dataf] = readhtk([datapath, blockname, filesep, whichdata, filesep, 'Ch', num2str(allchannels(j)), '.htk']);
    end
    
    %% Apply filters
    
    if highpass
        disp('...applying high pass filter');
        data_orig = data;
        [b,a] = butter(1,highpassfreq/(dataf/2),'high');
        data = filtfilt(b,a,data')';
        cleanvars.filters.highpassfreq = highpassfreq;
    end
    
    %% Remove bad channels
    
    if ~isempty(badchannels)
        disp('...removing bad channels');
        channels = 1:size(data,1);
        channels(badchannels) = [];
        data_cleanchannels = data;
        data_cleanchannels(badchannels,:) = [];
        
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
        % Save variables
        cleanvars.amplifiers.byorigchannel = amplifiersbyorigchannel;
        cleanvars.amplifiers.bythisout = amplifiers;
        
        
        %         data_bestchannels = data;
        %         data_bestchannels(cat(2,badchannels,iffychannels),:) = [];
    else
        data_cleanchannels = data;
    end
    
    
    %% referencing
    disp('...referencing data');
    
    switch reftype
        case 'all' % all data together, no grouping
            switch refsubtype
                case 'pca' % PCA with the proper number of components
                    disp(['...using ' num2str(bestPC) ' PCA components of all channels']);
                    datastd = mapstd(data_cleanchannels);
                    [u,~] = svd(datastd*datastd');
                    data_referenced = datastd - u(:,1:bestPC)*(datastd'*u(:,1:bestPC))';
                case 'avg'
                    disp('...using average of all channels');
                    data_referenced = data_cleanchannels - repmat(mean(data_cleanchannels,1),[size(data_cleanchannels,1) 1]);
                case 'med'
                    disp('...using median of all channels');
                    data_referenced = data_cleanchannels - repmat(median(data_cleanchannels,1),[size(data_cleanchannels,1) 1]);
            end
            
        case 'amp' % amplifiers are specified
            switch refsubtype
                case 'pca'
                    disp(['...using ' num2str(bestPC) ' PCA components of each amplifier']);
                    datastd = mapstd(data_cleanchannels);
                    data_referenced = zeros(size(data_cleanchannels,1), size(data_cleanchannels,2));
                    for j = 1:length(amplifiers)
                        datalocal = datastd(amplifiers{j},:);
                        [u,~] = svd(datalocal*datalocal');
                        data_referenced(amplifiers{j},:) = datalocal - u(:,1:bestPC) * (datalocal'*u(:,1:bestPC))';
                    end
                case 'avg'
                    disp('...using average of each amplifier');
                    data_referenced = zeros(size(data_cleanchannels));
                    for j = 1:length(amplifiers)
                        datalocal = data_cleanchannels(amplifiers{j},:);
                        data_referenced(amplifiers{j},:) = datalocal - repmat(mean(datalocal,1),[size(datalocal,1) 1]);
                    end
                case 'med'
                    disp('...using median of each amplifier');
                    data_referenced = zeros(size(data_cleanchannels));
                    for j = 1:length(amplifiers)
                        datalocal = data_cleanchannels(amplifiers{j},:);
                        data_referenced(amplifiers{j},:) = datalocal - repmat(median(datalocal,1),[size(datalocal,1) 1]);
                    end
            end
            
        case 'loc' % local
            disp('...starting local reference');
            warning('need to expand the local referencing to be possible on each amplifier');
            numchannels = size(data_cleanchannels,1);
            switch refsubtype
                case 'pca'
                    error('PCA for local reference is not currently supported');
                case 'avg'
                    disp('...using local average');
                    data_referenced = zeros(size(data_cleanchannels));
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
                        data_referenced(j,:) = thisdata - mean(datalocal,1);
                    end
                case 'med'
                    disp('...using local median');
                    data_referenced = zeros(size(data_cleanchannels));
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
                        data_referenced(j,:) = thisdata - mean(datalocal,1);
                    end
            end
        otherwise
            disp('...no reference selected, leaving data as is');
            data_referenced = data_cleanchannels;
    end
    
    
    
    
    %% epochs
    
    if ~isempty(badepochs)
        disp('...removing bad epochs');
        warning('needs to be fixed so that badepochs are considered by time, not by sample');
        data_cleanepochs = data_referenced;
        data_cleanepochs(:,badepochs) = [];
        
        warning('need to load and clip audio as well');
        %         aud_cleanepochs = aud;
        %         aud_cleanepochs(:,badepochs) = [];
        %
        data_cleanepochs = data_cleanepochs;
    else
        data_cleanepochs = data_referenced;
    end
    
    %% normalizing
    
    if normdata
        disp('...z-scoring data');
        data_cleaned = zscore(data_cleanepochs,[],2);
    else
        data_cleaned = data_cleanepochs;
    end
    
    
    %% Show the results if desired
    if visualizeresults
        figure;
        if highpass && normdata
            subplot(611);
            imagesc(data_orig); title('original');
            subplot(612)
            imagesc(data); title(['high pass filter at ' num2str(highpassfreq) ' Hz']);
            subplot(613)
            imagesc(data_cleanchannels); title('without bad channels');
            subplot(614)
            imagesc(data_referenced); title(['after ' fullreftype ' referencing']);
            subplot(615)
            imagesc(data_cleanepochs); title('without bad epochs');
            subplot(616)
            imagesc(data_cleaned); title('z-scored');
            
        elseif highpass
            subplot(511);
            imagesc(data_orig); title('original');
            subplot(512)
            imagesc(data); title(['high pass filter at ' num2str(highpassfreq) ' Hz']);
            subplot(513)
            imagesc(data_cleanchannels); title('without bad channels');
            subplot(514)
            imagesc(data_referenced); title(['after ' fullreftype ' referencing']);
            subplot(515)
            imagesc(data_cleaned); title('without bad epochs');
            
        elseif normdata
            subplot(511);
            imagesc(data); title('original');
            subplot(512)
            imagesc(data_cleanchannels); title('without bad channels');
            subplot(513)
            imagesc(data_referenced); title(['after ' fullreftype ' referencing']);
            subplot(514)
            imagesc(data_cleanepochs); title('without bad epochs');
            subplot(515)
            imagesc(data_cleaned); title('z-scored');
            
        else
            subplot(411);
            imagesc(data); title('original');
            subplot(412)
            imagesc(data_cleanchannels); title('without bad channels');
            subplot(413)
            imagesc(data_referenced); title(['after ' fullreftype ' referencing']);
            subplot(414)
            imagesc(data_cleaned); title('without bad epochs');
            
        end
        suptitle('Results of Cleaning Data');
        colormap jet
    end
    
    
    
    %% Now write the data
    
    CMRtag = 'cleaned';
    channels = allchannels;
    channels(badchannels) = [];
    
    % Check whether this block folder already exists; if so, double-check user wants to continue and overwrite
    if exist([datapath, blockname, filesep, dir_tag, CMRtag],'dir')
        processblock = input(['This block''s ' CMRtag ' folder already exists. Overwrite it? '],'s');
        if ~strcmpi(processblock,{'y', 'yes', '1'})
            disp('To avoid overwrite, revise block list and rerun');
            break
        end
    end
    
    % Make CMR directory
    disp(['...making ' dir_tag CMRtag ' directory']);
    mkdir([datapath, blockname, filesep, dir_tag, CMRtag]);
    
    % Write htk files for each specified data channel
    disp('...writing to file')
    for j = 1:length(channels) % loop to make htk for each channel
        writehtk([datapath, blockname, filesep, dir_tag, CMRtag, filesep, 'Ch', num2str(channels(j)), '.htk'], data_cleaned(j,:), dataf);
    end
    
    % Save the cleaning variables for info later
    disp('...saving cleaning variables');
    save([datapath, blockname, filesep CMRtag filesep 'cleaningvars.mat'],'cleanvars');
    
    disp(['Completed ' blockname]);
    
end

fprintf('Completed HTKtoCleanData \n');

end

