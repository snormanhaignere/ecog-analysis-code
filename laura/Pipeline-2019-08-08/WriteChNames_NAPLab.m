function [ chnames ] = WriteChNames_NAPLab(dpath,blocks,trodenames,trodelengths,loadchnames,startingpoint,trodeinds_backwards)
% [ chnames ] = WriteChNames_NAPLab(dpath,blocks,trodenames,trodelens,loadchnames,startingpoint,trodeinds_backwards)
%       Writes custom channel names; only necessary if channel names are not available in original data or there are errrors
%
%   Inputs:
%       datapath (optional): subject folder
%           datapath must contain raw data in an 'original' subdirectory
%               eg: datapath / original <- contains raw ECoG
%           default: current directory
%       blocks: specify which blocks to save these chnames to
%           format: array of integers,  eg: [1:4] 
%           default: all blocks in datapath / original
%       trodenames: cell array of channel name prefixes
%           format: cell array of strings
%       trodelengths: array of how many channels have each prefix
%           1:1 correspondence between trodenames and trodelengths
%           format: array of integers; 
%       loadchnames: if true, load existing channel names as a reference &
%           activate "keyboard" to allow user to investigate them
%           format: logical
%           default: false
%       startingpoint: which index of chnames to begin with
%           format: integer
%           default: 1 (beginning)
%
%       trodeinds_backwards: if true, indices will be applied backwards
%           (included due to technician errors in particular patients)
%           format: logical
%           default: false
%
%   Output:
%       chnames: cell array of resulting channel names

%% Params
if ~exist('loadchnames','var') || isempty(loadchnames)
    loadchnames = 0;
end
if ~exist('startingpoint','var') || isempty(startingpoint)
    startingpoint = 1;
end
if ~exist('trodeinds_backwards','var') || isempty(trodeinds_backwards)
    trodeinds_backwards = 0;
end

%% Load channel names if applicable
clear chnames;
if loadchnames
    load([dpath filesep 'original/B' num2str(blocks(1)) filesep 'allchnames_B' num2str(blocks(1)) '.mat'], 'chnames')
    disp('chnames already exists');
    keyboard;
else
    chnames = cell(1,numch);
end


%% Load trodenames into chnames
% Loop over trodenames and load into chnames
for i = 1:length(trodenames);
    
    if length(trodelengths) == 1
        trodelengths = repmat(trodelengths,1,length(trodenames));
    end
    thesenames = cell(1,trodelengths(i));
    for j = 1:length(thesenames);
        if ~trodeinds_backwards
            thesenames(j) = {[trodenames{i} num2str(j) ]}; % if proper order, just use the counter
        else
            thesenames(j) = {[trodenames{i} num2str(trodelengths(i) - j + 1)]}; % if backwards, fix
        end 
    end
    endpoint = startingpoint+trodelengths(i)-1; % the end of this channel set
    chnames(startingpoint: endpoint) = thesenames; % load the current names into chnames
    startingpoint = endpoint+1; % update the endpoint
    
end

disp(cat(1,chnames,num2cell(1:length(chnames)))');


%% Save 
savechnames = input('Save chnames? ','s');

for i = 1:length(blocks)
    if ~strcmpi(savechnames,{'y', 'yes', '1'})
    else
        savename = [dpath '/original/B' num2str(blocks(i)) filesep 'allchnames_B' num2str(blocks(i)) '.mat'];
        save(savename, 'chnames')
    end
end

end

