function [ channelind ] = getchannelinds(datapath,datatype)
%   [channelind] = getchannelinds(datapath,datatype)
%   Part of ECoG Pipeline.
%   Reads files in datapath and returns array of channel #s.
%
%   Currently works for TDT, need to add more options
%
% Written by Laura 10/2015


filenames = getfilenames(datapath);

channelind = [];
switch lower(datatype)
    case 'tdt'
        for i=1:length(filenames)
            chfind = strfind(filenames{i},'ch');
            if ~isempty(chfind) % if the file is a block folder (begins with 'B')
                currblockind = char(filenames{i}); % turn this into a char
                channelind(i) =  str2num(currblockind(chfind+2:end-4)); % now take the number and store
                
            end
        end
        
    case 'htk'
        for i = 1:length(filenames)
            chfind = strfind(filenames{i},'Ch');
            if ~isempty(chfind)
                currblockind = char(filenames{i});
                channelind(i) = str2num(currblockind(3:end-4));
            else
                chfind = strfind(filenames{i},'a');
                if length(chfind) == 1 && strcmp(filenames{i}(end-2:end),'htk')
                    currblockind = char(filenames{i});
                    channelind(i) = str2num(currblockind(2:end-4));
                end
            end
            
        end
end


channelind = sort(channelind);
channelind = channelind(channelind ~= 0);


end

