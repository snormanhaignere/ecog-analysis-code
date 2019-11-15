function StimOrder = makeStimFolder_Laura(datapath,blocks,stimpath,task,trials,portion)
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
%       portion: cell array indicating which portion of the StimOrder was
%           played if the experiment was interrupted or broken into two blocks;
%           portion{i} corresponds to block{i}
%
% Result:
%       Creates 'Stimulus' directory in each processed block in pt file.
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


disp(sprintf('\nMakeStimFolder'));

% Generate blocknames as cell array
for i = 1:length(blocks)
    blocknames{i} = ['B' num2str(blocks(i))];
end

% Default trials to the first trial
if ~exist('trials','var') || isempty(trials)
    trials = ones(1,length(blocks));
end

% Determine whether a portion will be taken
if ~exist('portion','var') || isempty(portion)
    isportion = 0;
else
    isportion = 1;
end

% Loop over blocks to make stimulus folder for each one
for i = 1:length(blocknames)
    
    % Display info and clear StimOrder in case it was preloaded (or from last block)
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
    
    % Method for finding StimOrder will change depending on task
    switch lower(task)
        case 'emotionalspeech' % this will have a patient-specific order due to randomization
            disp('...finding PtOrder from .mat file');
            behavfiles = getfilenames([datapath filesep 'behavioral' filesep blocknames{i} filesep]); % load block's behav files
            ptorderfile = [datapath filesep 'behavioral' filesep blocknames{i} filesep behavfiles{end}]; % find the mat file
            load(ptorderfile);
            resp = sortstruct(resp,'trial'); % sort by trial
            PtOrder = [resp.audiofilenumber]; % extract PtOrder
            StimOrder = {resp.audiofile};
            
            % Take out portion if applicable
            if isportion
                resp = resp(portion{i});
                PtOrder = PtOrder(portion{i});
            end
            
            
        case 'segmentation' % this will have patient-specific order due to randomization; may also be Presentation logfiles
            
            % Find log and mat files
            behavpath = [datapath filesep 'behavioral' filesep blocknames{i}];
            matfiles = dir([behavpath filesep '*.mat']);
            logfiles = dir([behavpath filesep '*.log']);
            
            if ~isempty(matfiles)
                disp('...finding PtOrder from .mat file');
                ptorderfile = [behavpath filesep matfiles(1).name];
                load(ptorderfile);
                StimOrder = audioplayed;
                
            elseif ~isempty(logfiles)
                
                disp('...finding StimOrder from LogFiles');
                
                % Find sort, and load logfile names
                logfiles = sortstruct(logfiles,'date'); % sort so that they're in chronological order
                logfilenames = {logfiles.name}; % get just the names
                
                % Load StimOrder
                StimOrder = {};
                for j = 1:length(logfilenames)
                    logfile = [behavpath filesep logfilenames{j}];
                    
                    % Read the log file
                    disp(['...reading from logfile ' logfilenames{j}]);
                    [~, out2] = importPresentationLog(logfile); %
                    
                    % Find the sounds and get the audiofile names for these
                    issound = strcmpi(out2.event_type,'Sound');
                    theseStim = out2.audiofile_str(issound)';
                    
                    % Add to the StimOrder
                    StimOrder = cat(2,StimOrder,theseStim);
                    
                end
            end
            
        case {'scrambling' 'pitch' 'naturalsounds' 'nyulocalizer', 'speech-tci-v2'}
            
            % Find log and mat files
            behavpath = [datapath filesep 'behavioral' filesep blocknames{i}];
            matfiles = dir([behavpath filesep '*.mat']);
            
            if ~isempty(matfiles)
                disp('...finding PtOrder from .mat file');
                StimOrder = {};
                for matfile = 1:length(matfiles)
                    ptorderfile = [behavpath filesep matfiles(matfile).name];
                    a = load(ptorderfile);
                    switch lower(task)
                        case {'naturalsounds'}
                            StimOrder = cat(2,StimOrder,a.allstim_for_this_run);
                        otherwise
                            try
                                StimOrder = cat(2,StimOrder,a.StimOrder(:)');
                            catch
                                keyboard
                            end
                    end
                    clear a;
                end
            end
            
            
        otherwise % For any other task, find and load appropriate stimorder by task/trial
            disp('...finding StimOrder');
            stimorderpath = [stimpath filesep task];
            stimorderfile = [stimorderpath filesep 'B' num2str(trials(i)) '_stimorder.mat'];
            load(stimorderfile);
            
    end
    
    % Take out portion if applicable
    if isportion
        StimOrder = StimOrder(portion{i});
    end
    
    % Generate Stimulus folder, write StimOrder to file inside
    disp('...writing to Stimulus folder');
    mkdir([datapath filesep 'processed' filesep blocknames{i} filesep, 'Stimulus']);
    save([datapath filesep 'processed' filesep blocknames{i} filesep, 'Stimulus' filesep 'StimOrder.mat'],'StimOrder');
    
    % Save PtOrder and sortresp if applicable
    if exist('PtOrder','var') && ~isempty(PtOrder)
        save([datapath filesep 'processed' filesep blocknames{i} filesep 'Stimulus' filesep 'PtOrder.mat'],'PtOrder');
        save([datapath filesep 'processed' filesep blocknames{i} filesep 'Stimulus' filesep 'resp.mat'],'resp');
    end
    
    disp(['Completed ' blocknames{i}]);
end

fprintf('\nCompleted MakeStimulusFolder \n');

end
