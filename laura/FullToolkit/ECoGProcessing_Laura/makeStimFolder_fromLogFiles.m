function makeStimFolder_fromLogFiles(datapath,blocks,stimpath,task)


fprintf('\nMakeStimFolder_fromLogFiles');

% Generate blocknames as cell array
blocknames = cell(1,length(blocks));
for i = 1:length(blocks)
    blocknames{i} = ['B' num2str(blocks(i))];
end

% if ~exist('portion','var') || isempty(portion)
%     isportion = 0;
% else
%     isportion = 1;
% end

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
    

    % If stimorder changes by patient, find ptorder
    if strcmpi(task,'Segmentation') || strcmpi(task,'SegmentationControl')
        disp('...finding StimOrder from LogFiles');
        
        % Find sort, and load logfile names
        logpath = [datapath filesep 'behavioral' filesep blocknames{i}];
        logfiles = dir([logpath filesep '*.log']); % load block's behav files
        logfiles = sortstruct(logfiles,'date'); % sort so that they're in chronological order
        logfilenames = {logfiles.name}; % get just the names
        
        StimOrder = {};
        for j = 1:length(logfilenames)
            logfile = [logpath filesep logfilenames{j}];
            
            % Read the log file
            disp(['...reading from logfile ' logfilenames{j}]);
            [~, out2] = importPresentationLog(logfile); % 
            
            % Find the sounds and get the audiofile names for these
            issound = strcmpi(out2.event_type,'Sound');
            theseStim = out2.audiofile_str(issound)';
            
            % Add to the StimOrder
            StimOrder = cat(2,StimOrder,theseStim);
            
        end
     
%         % Take out portion if applicable
%         if isportion 
%             resp = resp(portion{i});
%             PtOrder = PtOrder(portion{i});out2
%         end


    else
        error('Wrong function- this works only for segmentation at this time.');
    end
    
    
        
    % Generate Stimulus folder, write StimOrder to file inside
    disp('...writing to Stimulus folder');
    mkdir([datapath filesep 'processed' filesep blocknames{i} filesep, 'Stimulus']);
    save([datapath filesep 'processed' filesep blocknames{i} filesep, 'Stimulus' filesep 'StimOrder.mat'],'StimOrder');
    
    % Save logfiles if applicable so we can figure where it came from
    if exist('logfiles','var') && ~isempty(logfiles)
        save([datapath filesep 'processed' filesep blocknames{i} filesep 'Stimulus' filesep 'logfiles.mat'],'logfiles');
    end
    
    disp(['Completed ' blocknames{i}]);
end

fprintf('\nCompleted MakeStimulusFolder \n');


end
