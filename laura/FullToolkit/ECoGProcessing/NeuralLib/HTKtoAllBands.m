function HTKtoAllBands(basepath,whichdata,blocks,channels,dataf)
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
% Last Updated: 10/23/2015 by Laura; added option to use non-htkraw
% folders, whichdata input, and whichdata tag to created filenames
%
%% Set inputs, find path

fprintf('\nHTKtoAllBands\n');

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
    channels = getchannelinds([datapath 'B' num2str(blocks(1)) filesep 'htkraw'],'htk');
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
    disp(['...making ' dirtag 'all band directories']);
    mkdir([datapath, blockname, filesep, dirtag, 'delta']);
    mkdir([datapath, blockname, filesep, dirtag, 'theta']);
    mkdir([datapath, blockname, filesep, dirtag, 'alpha']);
    mkdir([datapath, blockname, filesep, dirtag, 'beta']);
    mkdir([datapath, blockname, filesep, dirtag, 'gamma']);
    mkdir([datapath, blockname, filesep, dirtag, 'highgamma']);
    
    
    % Write htk files for each specified data channel
    disp('...extracting all bands and writing to file')
    for j = channels % loop to make htk for each channel
        data = []; data_allbands = [];
        [data, oldf] = readhtk([datapath, blockname, filesep, whichdata, filesep, 'Ch', int2str(j), '.htk']);
        data_allbands = EcogExtractAllBands(data, oldf, dataf);
        writehtk([datapath, blockname, filesep, dirtag, 'delta' filesep, 'Ch', int2str(j), '.htk'], data_allbands(1,:), dataf);
        writehtk([datapath, blockname, filesep, dirtag, 'theta' filesep, 'Ch', int2str(j), '.htk'], data_allbands(2,:), dataf);
        writehtk([datapath, blockname, filesep, dirtag, 'alpha' filesep, 'Ch', int2str(j), '.htk'], data_allbands(3,:), dataf);
        writehtk([datapath, blockname, filesep, dirtag, 'beta' filesep, 'Ch', int2str(j), '.htk'], data_allbands(4,:), dataf);
        writehtk([datapath, blockname, filesep, dirtag, 'gamma' filesep, 'Ch', int2str(j), '.htk'], data_allbands(5,:), dataf);
        writehtk([datapath, blockname, filesep, dirtag, 'highgamma' filesep, 'Ch', int2str(j), '.htk'], data_allbands(6,:), dataf);
    end
    
    disp(['Completed ' blockname]);
    
end

fprintf('Completed HTKtoAllBands \n');

end

