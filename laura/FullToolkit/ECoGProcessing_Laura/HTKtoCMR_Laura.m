function CMRtag = HTKtoCMR(basepath,whichdata,blocks,channels,dataf,AvgorPCAorMed,amplifiers,badchannels)
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

fprintf('\nHTKtoCMR\n');

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
if ~exist('channels','var') || isempty(channels) % channels default to all
    channels = getchannelinds([datapath 'B' num2str(blocks(1)) filesep 'htkraw'],'htk');
end

if ~exist('PCAorAvg','var') || isempty(AvgorPCAorMed)
    AvgorPCAorMed = 'Avg';
end

if ~exist('amplifiers','var') || isempty(amplifiers)
    ampref = 0;
else
    ampref = 1;
end

if ~exist('badchannels','var') || isempty(badchannels)
    badchannels = [];
end


%% Loop over blocks, creating CMR HTK files

for i = blocks
    
    % Generate block name
    blockname = ['B' int2str(i)]; % create block string for use in files
    disp(['Processing ' blockname]);
    
    % Throw error if whichdata file isn't already there
    if ~exist([datapath blockname, filesep, whichdata],'dir')
        error(['Cannot find ' whichdata ' directory; check processed folder']);
    end
    
    %% Take Common Reference Average
    
    % Load data
    disp('...loading data');
    data = [];
    for j = 1:length(channels)
        [data(j,:),dataf] = readhtk([datapath, blockname, filesep, whichdata, filesep, 'Ch', num2str(channels(j)), '.htk']);
    end
    dataref = zeros(size(data));
    
    % Take Common Reference according to specified conditions
    disp('...finding CMR');
    switch AvgorPCAorMed
         
        case 'PCA' % use principle component analysis
            
            if ampref % reference by each amplifier
                
                datastd = mapstd(data);
                for j = 1:length(amplifiers)
                    dataamp = datastd(amplifiers{j},:);
                    [u,s] = svd(dataamp*dataamp');
                    dataref(amplifiers{j},:) = dataamp - u(:,1) * (dataamp'*u(:,1))';
                end
                
                CMRtag = 'CMRPCAamp';
                
            else % reference all together
                datastd = mapstd(data);
                [u,s] = svd(datastd*datastd');
                dataref = datastd - u(:,1:4)*(datastd'*u(:,1:4))';
                
                CMRtag = 'CMRPCA';
            end
            
        case 'Avg'
            
            if ampref
                
                for j = 1:length(amplifiers)
                    dataamp = data(amplifiers{j},:);
                    dataref(amplifiers{j},:) = dataamp - repmat(mean(dataamp,1),[size(dataamp,1) 1]);
                end
                
                CMRtag = 'CMRamp';
            else
                dataref = data - repmat(mean(data,1),[size(data,1) 1]);
                CMRtag = 'CMR';
            end
            
        case 'Med'
            
            if ampref
                
                for j = 1:length(amplifiers)
                    dataamp = data(amplifiers{j},:);
                    dataref(amplifiers{j},:) = dataamp - repmat(median(dataamp,1),[size(dataamp,1) 1]);
                end
                
                CMRtag = 'CMRMedamp';
            else
                dataref = data - repmat(median(data,1),[size(data,1) 1]);
                CMRtag = 'CMRMed';
            end
            
    end
       
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
        writehtk([datapath, blockname, filesep, dir_tag, CMRtag, filesep, 'Ch', num2str(channels(j)), '.htk'], dataref(j,:), dataf);
    end
    
    disp(['Completed ' blockname]);
    
end

fprintf('Completed HTKtoCMR \n');

end

