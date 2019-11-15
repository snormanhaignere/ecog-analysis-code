function evnt = NeuralFindEvent (DataPath, SoundPath, Subject, Blocks)
% findSyncEegGeneral function is used to locate the synchronization
% example: evnt = findSyncEegSingle('./B02/htkFiles/AN.htk','EEG_fs1_Single','~/Documents/MATLAB/nima_lab_local/GUI/CUSounds/','minda','B02',stim);
% points of each sentence.
% Input:
%   DataPath: a path to the "CUEEGxx" folder
%   SoundPath: the path to the sound files
%   Subject: name and ID of the subject for example "CUEEG1"
%   Blocks: cell array, containg Blocks names, like {'B01','B05'}

if exist([DataPath 'processed'],'dir')
    DataPath = [Datapath 'processed/'];
end

evntIndex = 1;
evnt = struct;
for cnt1 = 1:length(Blocks)
    display(['Block ' int2str(cnt1) ' Processing ...']);
    StimOrderPath = [DataPath Blocks{cnt1} filesep 'Stimulus/StimOrder.mat'];
    tmp= load(StimOrderPath);
    StimOrder=tmp.StimOrder;
    StimPath = [DataPath Blocks{cnt1} '/analog/a1.htk'];
    [audioRecord, audioRecordFreq] = readhtk(StimPath); % audio recording
    audioRecordFreq=audioRecordFreq;
    audioRecord = audioRecord(:)';
    
    syncPosition = [];
    stop = size(audioRecord,2);
    lastSyncPoint = 0;
    lastLength = 0;
    goToNext = 1; % condition flag on confidence, if confidence is big enough then go to the next sentence
    soundIndex = 1;
    while soundIndex <= length(StimOrder)  % loop the audio file names
        if goToNext == 1
            % Check if file exists
            SoundFile = [SoundPath filesep StimOrder{soundIndex}];
            % Downsample
            if exist('wavread')
                [w,fs_audio] = wavread(SoundFile);
            else
                [w,fs_audio] = audioread(SoundFile);
            end
            w2 = resample(w,audioRecordFreq,fs_audio);  %downsample to 2400 Hz 
            goToNext = 0;
        end
        
        % Cross-correlation
        windowSize = ceil(length(w2)/audioRecordFreq)+0.5;
        audioPart = audioRecord(lastSyncPoint+lastLength+1:min(lastSyncPoint+lastLength+windowSize*audioRecordFreq, stop)); % sampling rate is 10K, find the xcorr in the next 20s window
        c = xcorr(w2, audioPart);
        syncPoint = abs(length(audioPart)-(find(abs(c)==max(abs(c)))));
        syncPosition = [syncPosition; syncPoint+lastSyncPoint+lastLength];
        lastSyncPoint = syncPosition(end);
        lastLength = floor(length(w2)/20);
        
        % Calculate confidence
        waveform = (w2-mean(w2)) / std(w2);
        tmp = audioRecord(syncPosition(end)+1:syncPosition(end)+length(waveform));
        recordWaveform = (tmp-mean(tmp)) / std(tmp);
        confidence = abs((recordWaveform*waveform) / (waveform'*waveform));
        
        if confidence > 0.5
            goToNext = 1;
            lastLength =length(w2);
        end
        
        % add fields to evnt struct
        if goToNext == 1
            evnt(evntIndex).name = StimOrder{soundIndex};
            evnt(evntIndex).confidence = confidence;
            evnt(evntIndex).syncPosition = round((syncPosition(end)));
            evnt(evntIndex).startTime = (syncPosition(end)) / audioRecordFreq;
            evnt(evntIndex).stopTime = (syncPosition(end)+length(w2)) / audioRecordFreq;
            evnt(evntIndex).subject = Subject;
            evnt(evntIndex).block = Blocks{cnt1};
            evnt(evntIndex).trial = soundIndex;
            evnt(evntIndex).DataPath = DataPath;
            evnt(evntIndex).StimPath = SoundPath;
            soundIndex = soundIndex + 1;
            evntIndex = evntIndex + 1;
        end
    end
end
% check the correctness of the sync detection
figure('color','w');
plot(waveform); hold on; plot(recordWaveform,'r');

% save the evnt struct
%save(['evnt_',Subject,Blocks,task, '.mat'], 'evnt')
end

