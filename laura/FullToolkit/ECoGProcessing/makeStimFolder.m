function StimOrder = makeStimFolder_Laura(datapath,blocks,stimpath,task,trials)
% function makeStimulusFolder_Laura
%       []= makeStimulusFolder_Laura(datapath,blocks,stimpath,task,trials)
%       Part of ECoG Data Pipeline.
%       Finds appropriate StimOrder and writes to Stimulus folder for use
%       in NeuralFindEvent_Laura
%
% Inputs: 
%       datapath (optional): subject folder 
%           datapath must contain data in a 'processed' subdirectory
%               eg: datapath / processed <- contains processed ECoG
%           default: current directory
%       blocks (optional): specify which blocks to process
%           default: all blocks in datapath / original
%       stimpath: location of all task stimorders
%           stimpath must contain folders named by task 
%       task: name of task folder within stimpath
%           stimpath/task contains stimorder files for that task 
%           with the following naming scheme: B#_stimorder.mat
%       trials: array of which blocks correspond to which task trials
%           trials(i) is the trial number for blocks(i)
%           loads B#_stimorder.mat, where # is trials(i)
%
% Result:       
%       Creates 'Stimulus' directory in each processed block in pt file.
%       Creates .htk files for each audio and ECoG channel within block file.
%       
%       Resulting file structure: 
%           StimOrder: datapath / processed / B# / Stimulus / StimOrder.mat
%
% Hard-Coded Assumptions:
%       1. User has a folder with all task stimorders, organized in
%       subdirectories by task. This folder is 'stimorder', and the
%       appropriate subdirectory is 'task'.
%       2. Within each task folder, there are folders named
%       'B#_stimorder.mat' for each trial (# is trial #). Within each
%       stimorder file is a cell array with the names of the .wav files
%       played in order during that task trial.
%       3. All blocks being processed were the same task. (Must be run
%       multiple times to process separate tasks.)
%
% Written: 11/2015 by Laura Long


fprintf('\nMakeStimFolder');

% Generate blocknames as cell array
for i = 1:length(blocks)
    blocknames{i} = ['B' num2str(blocks(i))];
end


% Loop over blocks to make stimulus folder for each one
for i = 1:length(blocknames)
    
    disp(['Processing ' blocknames{i} '...']);
    
    clear 'StimOrder';
    
    
    % Check whether this block folder already exists; if so, double-check user wants to continue and overwrite
    if exist([datapath filesep 'processed' filesep blocknames{i} filesep, 'Stimulus'],'dir')
        processblock = input('This block''s Stimulus folder already exists. Overwrite it? ','s');
        if ~strcmpi(processblock,{'y', 'yes', '1'})
            disp('To avoid overwrite, revise block list and rerun');
            break
        end
    end
    
    % Find appropriate stimorder by task/trial and load
    disp('...finding StimOrder');
    stimorderpath = [stimpath filesep task];
    stimorderfile = [stimorderpath filesep 'B' num2str(trials(i)) '_stimorder.mat'];
    load(stimorderfile);
    
    
    % Generate Stimulus folder, write StimOrder to file inside
    disp('...writing to Stimulus folder');
    mkdir([datapath filesep 'processed' filesep blocknames{i} filesep, 'Stimulus']);
    save([datapath filesep 'processed' filesep blocknames{i} filesep, 'Stimulus' filesep 'StimOrder.mat'],'StimOrder');
    
    disp(['Completed ' blocknames{i}]);
end

fprintf('Completed MakeStimulusFolder \n');

end
