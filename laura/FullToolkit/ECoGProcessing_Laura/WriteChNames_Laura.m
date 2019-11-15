function [ chnames ] = WriteChNames_Laura(dpath,blocks,trodenames,trodelens,loadchnames,startingpoint,trodeinds_backwards)

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
    
    if length(trodelens) == 1
        trodelens = repmat(trodelens,1,length(trodenames));
    end
    thesenames = cell(1,trodelens(i));
    for j = 1:length(thesenames);
        if ~trodeinds_backwards
            thesenames(j) = {[trodenames{i} num2str(j) ]}; % if proper order, just use the counter
        else
            thesenames(j) = {[trodenames{i} num2str(trodelens(i) - j + 1)]}; % if backwards, fix
        end 
    end
    endpoint = startingpoint+trodelens(i)-1; % the end of this channel set
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

