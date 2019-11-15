function [filenames] = getfilenames(datapath)
%   [filenames] = getfilenames(datapath)
%   Returns cell array containing all filenames in the datapath
%
% Written by Laura 10/2015



files= dir(datapath); % Read all files in the folder

% Extract file names
filenames=cell(length(files),1); % preallocate filenames
for i=1:length(files)
    filenames{i}=files(i).name; % store filename
end

end