function evnt = NeuralFindEvent_Trigger_NAPLab(datapath, soundpath, pulsepath, subject, blocks, audchannel, troubleshoot)
% evnt = NeuralFindEvent_Trigger_NAPLab(datapath, soundpath, pulsepath, subject, blocks, audchannel, troubleshoot)
%       Creates evnt structure aligning .wav files to experimental data where only triggers are available
%
% Inputs: 
%       datapath (optional): subject folder (string) 
%           datapath must contain 'Stimulus' folder containing StimOrder.mat
%           within its 'processed' subdirectory
%               eg: datapath / Processed / B# / Stimulus
%           default: current directory
%       soundpath: folder of task .wav files (string) 
%           soundpath must contain all .wav files used in the experiment
%       pulsepath: .mat file containing original pulse path 
%           required if trigger.htk does not exist
%       subject: subject ID (will be used to name evnt, string)
%       blocks (optional): specify which blocks to process
%           format: array of integers OR cell array of block folder names
%               eg: [1:4] OR {'B1' 'B2' 'B3' 'B4'}
%           default: all blocks in datapath / original
%       audchannel (optional): specify which audio channel to use for
%           alignment (should be the computer input to TDT)
%           default: 1
%       troubleshoot: indices of stimuli for troubleshoot plotting
%           will pause, plot, and play alignment at each index
%           (if cell array, troublshoot{i} applies to block(i))
%
% Result:       
%       For each block, reads .wav files according to StimOrder
%       Finds correlation between .wav and experimentally-recorded sound
%       Finds start time of each .wav 
%       Displays sync for user approval
%       Asks user whether to save evnt file
%       
%       Resulting file: 
%           Evnt: datapath / processed / B# / evnt_subject_B#.mat
%
% Hard-Coded Assumptions:
%       1. StimOrder already exists in each block folder, with only the
%       stimuli actually presented (will not run for a partial block).
%       2. Triggers have already been processed, and there is only ONE
%       positive value for each stimulus. 


% This version adapted from Bahar's code by Laura 11/2015;
% - adjusted for new ECoG Data organization (processed folder)
% - cut second correlation iterance; it just screwed things up!
% - added option to save evnt if desired
% - removed rounding step in line 81; screwed up ECoG data because it could
% cut a file down to 0 samples
% - added lots of comments
% - revised blocks input to be valid as matrix or cell
% - evnt saves to each block folder
% - changed variable names (to lowercase, switched input order to be
% consistent with Convert/Out functions)
% 1/5/15
% - Changed default confidence to .6. 
% - Added option to select audio input; otherwise, looks for the first channel. 
% - Added detection of the audio tdt names; if it's not Aud_ or Aud1 returns an error. 
% - Added disp of current sound and confidence.

% Last Updated: 06/10/2019 by Laura
% - Cleaned up for NAPLab code library
% - Clarified addition of troubleshoot
% - Removed savespecgram option

%% Deal with inputs; set up evnt/evntIndex

fprintf('\nNeuralFindEvent\n');

% Check if datapath was entered; if not, default to current directory
if ~exist('datapath','var') || isempty(datapath)
    datapath = pwd; % if not specified, default to the current directory
end

% Blocks defaults to all blocks found in the original folder; populate blocknames
if ~exist('blocks','var') || isempty(blocks) % blocks default to all
    blocks = getblockinds([datapath filesep 'Processed' filesep]); %getblockinds finds vector of block numbers from datapath/original/
end
if iscell(blocks)
    blocknames = blocks;
else
    for i = 1:length(blocks)
        blocknames{i} = ['B' num2str(blocks(i))];
    end
end

% Check if audchannel was entered; if not, default to current directory
if ~exist('audchannel','var') || isempty(audchannel)
    audchannel = 1;
end

% Default troubleshoot to 0
if ~exist('troubleshoot','var') || isempty(troubleshoot) || troubleshoot == 0
    troubleshoot = [];
end



%% Loop over blocks
for i = 1:length(blocknames)
    
    disp(['Processing ' blocknames{i} '...']);
    evnt = struct; % set up a struct for evnt
    evntIndex = 1; % index for evnt structure
        
    % Read StimOrder
    StimOrderPath = [datapath filesep 'processed' filesep blocknames{i} filesep 'Stimulus' filesep 'StimOrder.mat']; % define path to stimorder mat file
    tmp = load(StimOrderPath); % load stimorder
    StimOrder = tmp.StimOrder; % define stimorder as StimOrder cell array
    StimPath = [datapath filesep 'processed' filesep blocknames{i} filesep 'analog' filesep 'a' num2str(audchannel) '.htk']; % define path to audio htk file
    
    
    troubleshootlogical = false(1,length(StimOrder));
    if iscell(troubleshoot), troubleshoot(troubleshoot{i}) = true;
    else troubleshootlogical(troubleshoot) = true; end
    
    
    %% Load trigger! 
    trigpath = [datapath filesep 'processed' filesep blocknames{i} filesep 'trigger' filesep];
    trigfile = [trigpath 'trigger.htk'];
    if exist(trigfile,'file')
        [trigger,fs_trig] = readhtk(trigfile);
    else
        disp('...matching trigger to pulse pattern');
        
        [t1,fs1] = readhtk([trigpath 't1.htk']);
        a = load(pulsepath); 
        t2 = a.trigform; fs2 = a.fs;
        [trigger,fs_trig] = matchTrigger(t1,fs1,t2,fs2);
        
        check = input('Accept this trigger alignment? ');
        if ~strcmpi(check,{'y', 'yes', '1'})
        else
            disp('...writing trigger.htk to file');
            save([trigpath 'trigger.mat'],'trigger','fs_trig');
            writehtk([trigpath 'trigger.htk'],trigger,fs_trig);
        end
        
    end
    
    if length(StimOrder) ~= sum(trigger)
        error('StimOrder length and number of trigger do not match');
    end
    
    triginds = find(trigger);
    
    %% Read the recorded experimental sound file
    if exist(StimPath,'file')
        [audioRecord, audioRecordFreq] = readhtk(StimPath); % reads the audio waveform and fs from the htk file
        audioRecord = audioRecord(:)'; % ensure the audio is a row vector
        [upsample,downsample] = rat(fs_trig/audioRecordFreq,1e-10);
        audioRecord = resample(audioRecord,upsample,downsample);
        wholewaveform = zeros(size(audioRecord));
    end
    
    % Set up sync variables
    syncPosition = []; % syncPosition is empty
    goToNext = 1; % condition flag for confidence; if confidence is big enough, go to the next sentence
    soundIndex = 1; % .wav file index (which stimulus in stimorder)
    
    %% Compare each .wav file in StimOrder with recorded sound
    while soundIndex <= length(StimOrder)  % loop the audio file names
        
        % Load the next .wav file if ready
        if goToNext == 1
            
            % Define the new sound file address
            SoundFile = [soundpath filesep StimOrder{soundIndex}];
            
            % If the soundfile is missing the .wav extension, add it
            if strcmpi(SoundFile(end-2:end),'wav') == 0 %
                SoundFile = [SoundFile, '.wav'];
            end
            
            % Read the audio file (with whichever function still exists for this computer)
            if exist('audioread')
                [w,fs_audio] = audioread(SoundFile);
            else
                [w,fs_audio] = wavread(SoundFile);
            end
            
            % Ensure it's one channel only
            if size(w,2) > 1
                w = w(:,1);
            end
            
            % Resample to match 
            [upsample,downsample] = rat(fs_trig/fs_audio,1e-10);
            w2 = resample(w,upsample,downsample);
            
            goToNext = 0; % set flag to 0
            
        end
        
        
        syncPosition = [syncPosition; triginds(soundIndex)];
        lastSyncPoint = syncPosition(end); % carry the chosen syncPoint forward

        
        %% Now find the correlation between the experimental sound and the .wav file!
        
        waveform = (w2-mean(w2)) / std(w2); % normalize .wav waveform
        
        if exist(StimPath,'file')
            % Calculate confidence
            tmp = audioRecord(syncPosition(end)+1:syncPosition(end)+length(waveform)); % find the recorded waveform
            recordWaveform = (tmp-mean(tmp)) / std(tmp); % normalize recorded waveform
            plot(recordWaveform); hold on; plot(waveform); hold off;
            confidence = abs((recordWaveform*waveform) / (waveform'*waveform)); % calculate correlation between normed .wav and recorded waveform
            
            % Troubleshooting will plot and play sounds
            if troubleshootlogical(soundIndex)
                disp(confidence)
                plot(waveform); hold on; plot(recordWaveform); hold off;
                clear sound;
                recordResamp = resample(recordWaveform,2,1);
                waveResamp = resample(waveform,2,1);
                soundsc(recordResamp,audioRecordFreq); soundsc(waveResamp,audioRecordFreq);
                keyboard;
            end
            
            % Go to next .wav
            goToNext = 1; % set flag to next
            disp(['Located sound ' num2str(soundIndex) ': Confidence ' num2str(confidence)]);
            lastLength =length(w2); % set last length to length of .wav
            wholewaveform(lastSyncPoint:lastSyncPoint+length(waveform)-1) = waveform;
            
        else
            confidence = NaN;
            goToNext = 1; % set flag to next
            disp(['Located sound ' num2str(soundIndex) ' by trigger']);
            lastLength =length(w2); % set last length to length of .wav
            wholewaveform(lastSyncPoint:lastSyncPoint+length(waveform)-1) = waveform;
        end
        
        %% Update the evnt structure if syncPosition was successfully found
        
        if goToNext == 1
            evnt(evntIndex).name = StimOrder{soundIndex}; % name of .wav
            evnt(evntIndex).confidence = confidence; % confidence of correlation
            evnt(evntIndex).syncPosition = round((syncPosition(end))); % sample that the audio begins
            evnt(evntIndex).startTime = (syncPosition(end)) / fs_trig; % start time of audio
            evnt(evntIndex).stopTime = (syncPosition(end)+length(w2)) / fs_trig; % stop time of audio
            evnt(evntIndex).subject = subject; % subject name
            evnt(evntIndex).block = blocknames{i}; % block number
            evnt(evntIndex).trial = soundIndex; % which sentence # (trial)
            evnt(evntIndex).DataPath = datapath; % where is the data?
            evnt(evntIndex).StimPath = soundpath; % where were the sounds?
            soundIndex = soundIndex + 1; % advance soundIndex
            evntIndex = evntIndex + 1; % advance evntIndex
        end
        
    end
    
    %% Finish block; display sync, ask about saving
    
    % Check the quality of sync by eye
    if exist(StimPath,'file')
        figure('color','w');
        wholerecord = (audioRecord-mean(audioRecord)) / std(audioRecord); % normalize recorded waveform
        plot(wholerecord,'b'); hold on; plot(wholewaveform,'r'); plot(trigger,'k','LineWidth',3);
        title('Record: Blue      From File: Red      Trigger: Black');
    else
        figure('color','w')
        plot(wholewaveform,'r'); hold on; plot(trigger,'k','LineWidth',2);
        title('From File: Red      Trigger: Black');
    end
    
    
    % Ask about saving
    saveevnt = input('Save evnt for this block? ','s');
    if ~strcmpi(saveevnt,{'y', 'yes', '1'})
    else
        save([datapath filesep 'processed' filesep blocknames{i} filesep 'evnt_' subject '_' blocknames{i}, '.mat'], 'evnt')
        disp(['Saved to ' datapath filesep 'processed' filesep blocknames{i} filesep 'evnt_' subject '_' blocknames{i}, '.mat']);
    end
    
    savewaveform = input('Save waveform for this block? ','s');
    if ~strcmpi(savewaveform,{'y', 'yes', '1'})
    else
        w = wholewaveform;
        savename = [datapath filesep 'processed' filesep blocknames{i} filesep 'Stimulus' filesep 'waveform.htk'];
        writehtk(savename,w,fs_trig);
        disp(['Saved to ' savename]);
    end
    
    % Notify the block is done
    disp(['Completed ' blocknames{i}]);
    
end

fprintf('Completed NeuralFindEvent \n');

end
