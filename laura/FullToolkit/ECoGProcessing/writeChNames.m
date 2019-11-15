%% Params
subject = '040_LIJ116'
loadchnames = 0; % decide whether to load already-written
numch = 156; % how many channels total
blocks = [1:3 11 14 17]; % which blocks do these channel apply to

%% Load channel names if applicable
clear chnames;
if loadchnames
    load(['/Users/LauraLong/Documents/Lab/ECoG Data' filesep subject filesep 'original/B' num2str(blocks(1)) filesep 'allchnames_B' num2str(blocks(1)) '.mat'], 'chnames')
    keyboard;
else
    chnames = cell(1,numch);
end


%% Define channel names

% list of trode names
trodename = {'LFo' 'out' 'LFa' 'out' 'LFm' 'out' 'LFi' 'out' 'LFp' 'out' 'LTx' 'out' 'LTs' 'LDa' ...
    'out' 'LDh' 'out' 'RFo' 'out' 'RFa' 'out' 'RFi' 'RTx' 'RTs' 'RTp' 'RDa' 'RDh'}; 
% how many channels each name applies to; if scalar, assumes it's the same for all
trodelen = [14 2 14 2 13 3 15 1 14 2 13 3 16 12 4 12 4 14 18 13 3 14 8 11 7 13 11]; 
startingpoint = 1; % where to start loading these new names into chnames
trodeinds_backwards = 0; % if the electrode indices are backwards (true for LIJ114)

% Loop over trodenames and load into chnames
for i = 1:length(trodename);
    
    
    if length(trodelen) == 1
        trodelen = repmat(trodelen,1,length(trodename));
    end
    
    thesenames = cell(1,trodelen(i));
    for j = 1:length(thesenames);
        
        if ~trodeinds_backwards
            thesenames(j) = {[trodename{i} num2str(j) ]}; % if proper order, just use the counter
        else
            thesenames(j) = {[trodename{i} num2str(trodelen(i) - j + 1)]}; % if backwards, fix
        end
        
    end
    
    endpoint = startingpoint+trodelen(i)-1; % the end of this channel set
    chnames(startingpoint: endpoint) = thesenames; % load the current names into chnames
    startingpoint = endpoint+1; % update the endpoint
    
end

disp(cat(1,chnames,num2cell(1:length(chnames)))');


%% Save
% Ask about saving
savechnames = input('Save chnames? ','s');

for i = 1:length(blocks)
    if ~strcmpi(savechnames,{'y', 'yes', '1'})
    else
        savename = ['/Users/LauraLong/Documents/Lab/ECoG Data/' subject '/original/B' num2str(blocks(i)) filesep 'allchnames_B' num2str(blocks(i)) '.mat'];
        save(savename, 'chnames')
    end
end