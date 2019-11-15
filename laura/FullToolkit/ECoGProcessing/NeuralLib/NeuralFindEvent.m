function evnt = NeuralFindEvent (DataPath, SoundPath, Subject, Blocks)
% findSyncEegGeneral function is used to locate the synchronization
% example: evnt = findSyncEegSingle('./B02/htkFiles/AN.htk','EEG_fs1_Single','~/Documents/MATLAB/nima_lab_local/GUI/CUSounds/','minda','B02',stim);
% points of each sentence.
% Input:
%   DataPath: a path to the "CUEEGxx" folder
%   SoundPath: the path to the sound files
%   Subject: name and ID of the subject for example "CUEEG1"
%   Blocks: cell array, containg Blocks names, like {'B01','B05'}
% updated by Bahar  June 2015


evntIndex = 1;
evnt = struct;
for cnt1 = 1:length(Blocks)
    tmplast=0;
    display(['Block ' int2str(cnt1) ' Processing ...']);
    StimOrderPath = [DataPath Blocks{cnt1} filesep 'Stimulus/StimOrder.mat'];
    tmp= load(StimOrderPath);
    StimOrder=tmp.StimOrder;
    StimPath = [DataPath Blocks{cnt1} '/analog/a1.htk'];
    [audioRecord, aduioRecordFreq] = readhtk(StimPath); % audio recording
    audioRecordFreq=aduioRecordFreq;
    audioRecord = audioRecord(:)';
    
    syncPosition = [];
    stop = size(audioRecord,2);
    lastSyncPoint = 0;
    lastLength = 0;
    goToNext = 1; % condition flag on confidence, if confidence is big enough then go to the next sentence
    soundIndex = 1;
    while soundIndex <=length(StimOrder)  % loop the audio file names
        if soundIndex == 145
            disp('');   
        end
        if goToNext == 1
            % Check if file exists
%             x = [SoundPath filesep StimOrder{soundIndex}];
            SoundFile = [SoundPath filesep StimOrder{soundIndex}];
            
            if strcmp(SoundFile(end-2:end),'wav') == 0
                SoundFile = [SoundFile, '.wav'];
            end
    
            
            % Downsample
            if exist('audioread')
                [w,fs_audio] = audioread(SoundFile);
            else
                [w,fs_audio] = wavread(SoundFile);
       
            end
            w3=resample(w,aduioRecordFreq,fs_audio);
            w=w(1:end-rem(length(w),aduioRecordFreq));
            w2 = resample(w,aduioRecordFreq,fs_audio);  %downsample to 2400 Hz 
            goToNext = 0;
        end
        
        % FROM OLD FUNCTION BELOW
        
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
        
        if confidence > 0.1
            goToNext = 1;
            lastLength =length(w2);
        end        
        
        % FROM OLD FUNCTION ABOVE
        
        
        % Cross-correlation
        windowSize = ceil(length(w2)/aduioRecordFreq)+0.5;
        audioPart = audioRecord(lastSyncPoint+lastLength+1:min(lastSyncPoint+lastLength+windowSize*aduioRecordFreq, stop)); % sampling rate is 10K, find the xcorr in the next 20s window
        c = xcorr(w2, audioPart);
        syncPoint = abs(length(audioPart)-(find(abs(c)==max(abs(c)))));
        syncPosition = [syncPosition; syncPoint+lastSyncPoint+lastLength];
        
        lastSyncPoint = syncPosition(end)-2;
        while  tmplast>=lastSyncPoint
            lastSyncPoint=lastSyncPoint+1;
        end
        
        tmplast=lastSyncPoint;
        lastLength = 1;%ceil(length(w2)/600); %% 60 is step size and it should change based on the length of the waveform
        
        % Calculate confidence
        waveform = (w2-mean(w2)) / std(w2);
        % TODO: Check here, mike
%         if((syncPosition(end)+length(waveform)) > length(audioRecord))
%             confidence = 0;
% %             disp('Zero Confidence');
%         else
%             tmp = audioRecord(syncPosition(end)+1:syncPosition(end)+length(waveform) );
% %         tmp = audioRecord(syncPosition(end)+1:min(syncPosition(end)+length(waveform),length(audioRecord)) );
%         recordWaveform = (tmp-mean(tmp)) / std(tmp);
%         confidence = abs((recordWaveform*waveform) / (waveform'*waveform));
%         if confidence > 0.5
%             goToNext = 1;
%             lastLength =length(w2);
%         end
%         end
        tmp = audioRecord(syncPosition(end)+1:syncPosition(end)+length(waveform) );
%         tmp = audioRecord(syncPosition(end)+1:min(syncPosition(end)+length(waveform),length(audioRecord)) );
        recordWaveform = (tmp-mean(tmp)) / std(tmp);
        len_w=round(length(recordWaveform));
        confidence = abs((recordWaveform(1:len_w)*waveform(1:len_w)) / (waveform(1:len_w)'*waveform(1:len_w)))
%         confidence2 = abs((recordWaveform(end-1000:end)*waveform(end-1000:end)) / (waveform(end-1000:end)'*waveform(end-1000:end)));
%         confidence=(confidence+confidence2)/2;
        if confidence > 0.25
            goToNext = 1;
            lastLength =length(w2);
        end
        
        
        % add fields to evnt struct
        if goToNext == 1
            evnt(evntIndex).name = StimOrder{soundIndex};
            evnt(evntIndex).confidence = confidence;
            evnt(evntIndex).syncPosition = round((syncPosition(end)));
            evnt(evntIndex).startTime = (syncPosition(end)) / aduioRecordFreq;
            evnt(evntIndex).stopTime = (syncPosition(end)+length(w3)) / aduioRecordFreq;
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
%save(['evnt_',Subject, '_', Blocks{1}, '.mat'], 'evnt')
end

