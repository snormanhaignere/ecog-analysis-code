function newtag = HTKtoFreqBand_NAPLab(basepath,whichdata,blocks,channels,dataf,freqRange,getenvelope,datatag)
% function HTKtoHighGamma
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
%       dataf (optional): desired output fs
%           default: 1000Hz
%       freqRange: range of frequencies to extract (e.g. [0 70])
%       getenvelope: logical of whether to take envelope of signal
%       datatag: name of data folder to save
% Output:
%       newtag: name of created folder
%
% Result:
%       Creates 'datatag' directory in each processed block file.
%       Extracts frequency range from files in basepath / processed / whichdata.
%       Writes a new .htk file with high gamma.
%
%       Resulting file structure:
%           From htkraw: basepath / processed / B# / datatag / Ch#.htk
%           Data: basepath / processed / B# / whichdata_datatag / Ch#.htk
%
% Hard-Coded Assumptions:
%       1. Default data channels are all channels in the first block
%       2. Original data is in htkraw folder.
%
%
% Written: 10/2015 by Laura
% Updates: 10/23/2015 by Laura; added option to use non-htkraw
% folders, whichdata input, and whichdata tag to created filenames
%
% Last Updated: 06/10/2019 by Laura; added newtag and updating comments for
%               NAPLab code package
%
%% Set inputs, find path

fprintf('\nHTKtoFreqBand\n');

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
    dirtag = []; % no tag if converting from raw
else
    dirtag = [whichdata '_']; % include tag if using already processed data
end
newtag = [dirtag datatag];

% Set parameter defaults

% Data sampling rate defaults to 1000 (aud is hard coded to 2400)
if ~exist('dataf','var') || isempty(dataf) % sampling frequency default to 1000
    dataf = 1000;
end

% Blocks defaults to all blocks found in the original folder
if ~exist('blocks','var') || isempty(blocks) % blocks default to all
    blocks = getblockinds(datapath); %getblockinds finds vector of block numbers from basepath/processed/
end

if ~exist('blocks','var') || isempty(blocks)
    blocks = [];
end
if iscell(blocks) % If blocks input is matrix, turn it into a cell array
    blocknames = blocks;
else
    for i = 1:length(blocks)
        blocknames{i} = ['B' num2str(blocks(i))];
    end
end

% Channels defaults to all channels found in the first block folder
if ~exist('channels','var') || isempty(channels) % channels default to all
    channels = getchannelinds([datapath 'B' num2str(blocks(1)) filesep whichdata],'htk');
end

%% Loop over blocks, creating high gamma HTK files

for i = blocks
    
    % Generate block name
    blockname = ['B' int2str(i)]; % create block string for use in files
    disp(['Processing ' blockname]);
    
    % Throw error if whichdata file isn't already there
    if ~exist([datapath blockname, filesep, whichdata],'dir')
        error(['Cannot find ' whichdata ' directory; check processed folder']);
    end
    
    % Check whether this block folder already exists; if so, double-check user wants to continue and overwrite
    if exist([datapath, blockname, filesep, dirtag, 'allbands'],'dir')
        processblock = input('This block''s all bands folder already exists. Overwrite it? ','s');
        if ~strcmpi(processblock,{'y', 'yes', '1'})
            disp('To avoid overwrite, revise block list and rerun');
            break
        end
    end
    
    % Make high gamma directory
    disp(['...making ' dirtag ' ' datatag ' directories']);
    mkdir([datapath, blockname, filesep, dirtag, datatag]);
    
    % Write htk files for each specified data channel
    disp(['...extracting ' datatag ' and writing to file'])
    for j = channels % loop to make htk for each channel
        data = []; data_newband = [];
        [data, oldf] = readhtk([datapath, blockname, filesep, whichdata, filesep, 'Ch', int2str(j), '.htk']);
        data_newband = EcogExtractFreqBand(data, oldf, dataf, freqRange, getenvelope);
        writehtk([datapath, blockname, filesep, dirtag, datatag, filesep, 'Ch', int2str(j), '.htk'], data_newband, dataf);
    end
    
    disp(['Completed ' blockname]);
    
end

fprintf(['Completed HTKtoFreqBand with tag ' datatag ' \n']);

end

