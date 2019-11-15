function evnt = NeuralFindEvent_Audio_NAPLab(datapath, soundpath, subject, blocks, audchannel, confidencethreshold, applyhilbert, checkmanually, troubleshoot, stepsize)
% evnt = NeuralFindEvent_Audio_NAPLab(datapath, soundpath, subject, blocks, audchannel, confidencethreshold, applyhilbert, checkmanually, troubleshoot, stepsize)
%       Creates evnt structure aligning .wav files to audio recorded during experiment.
%
% Inputs:
%       datapath (optional): subject folder (string)
%           datapath must contain 'Stimulus' folder containing StimOrder.mat
%           within its 'processed' subdirectory
%               eg: datapath / Processed / B# / Stimulus
%           default: current directory
%       soundpath: folder of task .wav files (string)
%           soundpath must contain all .wav files used in the experiment
%       subject: subject ID (will be used to name evnt, string)
%       blocks (optional): specify which blocks to process
%           format: array of integers OR cell array of block folder names
%               eg: [1:4] OR {'B1' 'B2' 'B3' 'B4'}
%           default: all blocks in datapath / original
%       audchannel: specify which audio channel to use for alignment
%           default: 1
%       confidencethreshold: correlation threshold for each stimulus
%           can either be a scalar, which will apply to all stimuli, 
%           or a vector with one value for each stimulus
%           (if cell array, confidencethreshold{i} applies to block(i))
%       applyhilbert: indices of stimuli that should be evaluated with
%           Hilbert transform (if failure to align due to phase issues)
%           (if cell array, confidencethreshold{i} applies to block(i))
%       checkmanually: indices of stimuli to adjust manually
%           (if cell array, confidencethreshold{i} applies to block(i))
%       troubleshoot: indices of stimuli for troubleshoot plotting
%           will pause, plot, and play alignment at each index
%           (if cell array, troublshoot{i} applies to block(i))
%       stepsize: integer of what fraction of stimulus length to step forward
%           default: 20
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
% - Clarified addition of Hilbert, manual check, troubleshoot,
% 	variable confidence threshold.
% - Added ability to change stepsize.


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
if ~exist('audchannel','var') || isempty(audchannel)
    audchannel = 1;
end
if ~exist('confidencethreshold','var') || isempty(confidencethreshold)
    confidencethreshold = 0.6;
end
if ~exist('needhilb','var') || isempty(applyhilbert) || applyhilbert == 0
    applyhilbert = [];
end
if ~exist('checkmanually','var') || isempty(checkmanually) || checkmanually == 0
    checkmanually = [];
end
if ~exist('troubleshoot','var') || isempty(troubleshoot) || troubleshoot == 0
    troubleshoot = [];
end
if ~exist('stepsize','var') || isempty(stepsize)
    stepsize = 20;
end


%% Loop over blocks
for i = 1:length(blocknames)
    
    disp(['Processing ' blocknames{i} '...']);
    evnt = struct; % set up a struct for evnt
    evntIndex = 1; % index for evnt structure
    
    % Read StimOrder from file
    StimOrderPath = [datapath filesep 'processed' filesep blocknames{i} filesep 'Stimulus' filesep 'StimOrder.mat']; % define path to stimorder mat file
    tmp = load(StimOrderPath); % load stimorder
    StimOrder = tmp.StimOrder; % define stimorder as StimOrder cell array
    StimPath = [datapath filesep 'processed' filesep blocknames{i} filesep 'analog' filesep 'a' num2str(audchannel) '.htk']; % define path to audio htk file
    
    % Set up vector of confidence thresholds
    if iscell(confidencethreshold)
        thisconfidencethreshold = confidencethreshold{i};
    else
        thisconfidencethreshold = confidencethreshold;
    end
    if length(thisconfidencethreshold) == 1
        thisconfidencethreshold = repmat(thisconfidencethreshold,1,length(StimOrder));
    elseif length(thisconfidencethreshold) ~= length(StimOrder)
        error('number of confidence thresholds does not match number of stimuli');
    end
    
    % Set up logicals for Hilbert transform and manual checking
    hilbertlogical = false(1,length(StimOrder));
    if iscell(applyhilbert), hilbertlogical(applyhilbert{i}) = true;
    else hilbertlogical(applyhilbert) = true; end
    manualchecklogical = false(1,length(StimOrder));
    if iscell(checkmanually), manualchecklogical(checkmanually{i}) = true;
    else manualchecklogical(checkmanually) = true; end
    troubleshootlogical = false(1,length(StimOrder));
    if iscell(troubleshoot), troubleshoot(troubleshoot{i}) = true;
    else troubleshootlogical(troubleshoot) = true; end
    
    % Read the recorded experimental sound file
    [audioRecord, audioRecordFreq] = readhtk(StimPath); % reads the audio waveform and fs from the htk file
    audioRecord = audioRecord(:)'; % ensure the audio is a row vector
    wholewaveform = zeros(size(audioRecord));
    
    % Set up sync variables
    syncPosition = []; % syncPosition is empty
    stop = size(audioRecord,2); % set stop var for once audioRecord is over
    lastSyncPoint = 0; % first syncPoint is 0
    lastLength = 0; % first length is 0
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
            try
                if exist('audioread')
                    [w,fs_audio] = audioread(SoundFile);
                else
                    [w,fs_audio] = wavread(SoundFile);
                end
            catch
                keyboard
            end
            
            % Resample to match
            [upsample,downsample] = rat(audioRecordFreq/fs_audio,1e-10);
            w2 = resample(w,upsample,downsample);
            w2 = w2(:,1);
            
            goToNext = 0; % set flag to 0
            
        end
        
        
        %% Now find the correlation between the experimental sound and the .wav file!
        
%         windowSize = ceil(length(w2)/audioRecordFreq)+0.5; % define comparison window size as the duration of the sound + .5 s
%         audioPart = audioRecord(lastSyncPoint+lastLength+1 : min(lastSyncPoint+lastLength+windowSize*audioRecordFreq, stop)); % clip part of the audio for comparison: between the last sync completion and the next windowSize (or the end of the file, whichever comes first)
        moveForward = max(lastLength-ceil(1e-3*audioRecordFreq),0); % Move forward to 1ms before the end of last segment
        moveForward = ceil(lastLength*0.99);
        windowSize = ceil(length(w2)/audioRecordFreq)+0.5; % define comparison window size as the duration of the sound + .5 s
        %audioPart = audioRecord(lastSyncPoint+lastLength+1 : floor(min(lastSyncPoint+lastLength+windowSize*audioRecordFreq, stop))); % clip part of the audio for comparison: between the last sync completion and the next windowSize (or the end of the file, whichever comes first)
        audioPart = audioRecord(lastSyncPoint+moveForward+1 : floor(min(lastSyncPoint+moveForward+windowSize*audioRecordFreq, stop)));

        
        
        % Find cross-correlation between recorded and actual
        c = xcorr(w2, audioPart);
        syncPoint = abs(length(audioPart)-(find(abs(c)==max(abs(c))))); % syncPoint in the window- where is the correlation best?
        syncPosition = [syncPosition; syncPoint+lastSyncPoint+moveForward]; % syncPosition in the entire file
        lastSyncPoint = syncPosition(end); % carry the chosen syncPoint forward
        lastLength = floor(length(w2)/stepsize); % step the length forward (to be used if it doesn't pass confidence threshold)
        
        % Calculate confidence
        waveform = (w2-mean(w2)) / std(w2); % normalize .wav waveform
        tmp = audioRecord(syncPosition(end)+1:syncPosition(end)+length(waveform)); % find the recorded waveform
        recordWaveform = (tmp-mean(tmp)) / std(tmp); % normalize recorded waveform
        
        if hilbertlogical(soundIndex)
            confidence = corr(abs(hilbert(recordWaveform))',abs(hilbert(waveform)));
            if troubleshootlogical(soundIndex)
                plot(waveform); hold on; plot(recordWaveform); hold off;
            end
        else
            confidence = abs((recordWaveform*waveform) / (waveform'*waveform)); % calculate correlation between normed .wav and recorded waveform
            if troubleshootlogical(soundIndex)
                plot(waveform); hold on; plot(recordWaveform); hold off;
            end
        end
        
        if troubleshootlogical(soundIndex)
            disp(confidence);
            plot(waveform); hold on; plot(recordWaveform); legend('from file','from recording');
            hold off;
        end
        
        
        % If the correlation was high enough, go to next .wav
        if confidence > thisconfidencethreshold(soundIndex) % Currently set to .6; can discuss to determine ideal threshold.
            
            % Adjust manually if necessary
            if manualchecklogical(soundIndex)

                readytomoveon = 0; firsttry = 1;
                while ~readytomoveon
                    
                    % Plot and select a better alignment location
                    figure;
                    plot(waveform); hold on; plot(recordWaveform); legend('from file','from recording'); title(confidence);
                    hold off;
                    if firsttry
                        needtoadjust = 'Adjust manually?'
                        if ~needtoadjust, break, end
                    end
                    [newSyncPoint,~] = ginput(1);
                    
                    % Update the syncpoint given this input
                    syncPosition(end) = ceil(syncPosition(end)+newSyncPoint);
                    lastSyncPoint = syncPosition(end); % carry the chosen syncPoint forward
                    lastLength = floor(length(w2)/stepsize); % step the length forward (to be used if the confidence isn't good enough)
                    
                    % Replot and calculate confidence
                    waveform = (w2-mean(w2)) / std(w2); % normalize .wav waveform
                    tmp = audioRecord(syncPosition(end)+1:syncPosition(end)+length(waveform)); % find the recorded waveform
                    recordWaveform = (tmp-mean(tmp)) / std(tmp); % normalize recorded waveform
                    if hilbertlogical(soundIndex)
                        confidence = corr(abs(hilbert(recordWaveform))',abs(hilbert(waveform)));
                    else
                        confidence = abs((recordWaveform*waveform) / (waveform'*waveform)); % calculate correlation between normed .wav and recorded waveform
                    end
                    figure;
                    plot(waveform); hold on; plot(recordWaveform); legend('from file','from recording'); title(confidence);
                    hold off;
                    
                    % Check if this was a good enough adjustment
                    readytomoveon = input('Accept and move on? ');
                end
                clear readytomoveon;
            end
            
            % In troubleshooting mode, play the sound and pause for the user to assess
            if troubleshoot
                clear sound;
                soundsc(recordWaveform,audioRecordFreq); soundsc(waveform,audioRecordFreq);
                plot(waveform); hold on; plot(recordWaveform);
                keyboard;
            end
            
            % Update flag and report result (warn if used Hilbert); update length and full waveform 
            goToNext = 1; % set flag to next
            disp(['Located sound ' num2str(soundIndex) ': Confidence ' num2str(confidence)]);
            if hilbertlogical(soundIndex)
                warning('Used Hilbert Transform');
            end
            lastLength =length(w2); % set last length to length of .wav
            wholewaveform(lastSyncPoint:lastSyncPoint+length(waveform)-1) = waveform;
            
        end
        
        
        %% Update the evnt structure if syncPosition was successfully found
        
        if goToNext == 1
            evnt(evntIndex).name = StimOrder{soundIndex}; % name of .wav
            evnt(evntIndex).confidence = confidence; % confidence of correlation
            evnt(evntIndex).syncPosition = round((syncPosition(end))); % sample that the audio begins
            evnt(evntIndex).startTime = (syncPosition(end)) / audioRecordFreq; % start time of audio
            evnt(evntIndex).stopTime = (syncPosition(end)+length(w2)) / audioRecordFreq; % stop time of audio
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
    figure('color','w');
    wholerecord = (audioRecord-mean(audioRecord)) / std(audioRecord); % normalize recorded waveform
    plot(wholerecord,'b'); hold on; plot(wholewaveform,'r');
    title('Record: Blue      From File: Red');
    
    % Notify the block is done
    disp(['Completed ' blocknames{i}]);
    
    % Ask about saving
    saveevnt = input('Save evnt for this block? ','s');
    if ~strcmpi(saveevnt,{'y', 'yes', '1'})
    else
        save([datapath filesep 'Processed' filesep blocknames{i} filesep 'evnt_' subject '_' blocknames{i}, '.mat'], 'evnt','thisconfidencethreshold')
        disp(['Saved to ' datapath filesep 'Processed' filesep blocknames{i} filesep 'evnt_' subject '_' blocknames{i}, '.mat']);
    end
    
    
end

fprintf('Completed NeuralFindEvent \n');

end



