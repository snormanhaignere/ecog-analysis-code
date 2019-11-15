function [blockind] = getblockinds(datapath)
%   [blockind] = getblockinds(datapath)
%    Part of ECoG Pipeline. 
%    Assumes datapath contains directories named 'B#' and returns array of block #s.
%
% Written by Laura 10/2015

filenames = getfilenames(datapath); % Get file names

% Loop through files; store block numbers in vector
blockind = [];
for i=1:length(filenames)
    if strfind(filenames{i},'B') == 1 % if the file is a block folder (begins with 'B')
        currblockind = char(filenames{i}); % turn this into a char
        blockind(i) =  str2num(currblockind(2:end)); % now take the number and store
        
    end
end

% Sort and remove any 0s
blockind = sort(blockind);
blockind = blockind(blockind ~= 0);

end