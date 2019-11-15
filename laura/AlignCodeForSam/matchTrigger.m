function [singletrigger,fs] = matchTrigger(trigRecord,trigRecordFreq,trigw,fs_trig,thresh,troubleshoot,manualcheck)
%

if ~exist('thresh','var') || isempty(thresh)
    thresh = .9;
end
if ~exist('troubleshoot','var') || isempty(troubleshoot)
    troubleshoot = 0;
end
if ~exist('manualcheck','var') || isempty(manualcheck)
    manualcheck = 0;
end

trigw = resample(trigw,trigRecordFreq,fs_trig);
fs = trigRecordFreq;

lastSyncPoint = 0;
lastLength = 0;
stop = size(trigRecord,2)-400;
syncPosition = [];
trigfoundlocs = [];

allconfidence = [];
windowSize = ceil(length(trigw)/trigRecordFreq)+0.5; % define comparison window size as the duration of the sound + .5 s

while lastSyncPoint+lastLength+windowSize*trigRecordFreq < stop
    
    audioPart = trigRecord(lastSyncPoint+lastLength+1 : min(lastSyncPoint+lastLength+windowSize*trigRecordFreq, stop)); % clip part of the audio for comparison: between the last sync completion and the next windowSize (or the end of the file, whichever comes first)
    
    % Find cross-correlation between recorded and actual
    c = xcorr(trigw, audioPart);
    syncPoint = abs(length(audioPart)-(find(abs(c)==max(abs(c))))); % syncPoint in the window- where is the correlation best?
    syncPosition = [syncPosition; syncPoint+lastSyncPoint+lastLength]; % syncPosition in the entire file
    lastSyncPoint = syncPosition(end); % carry the chosen syncPoint forward
    lastLength = floor(length(trigw)/20); % step the length forward 20 (to be used if the confidence isn't good enough)
    
    % Calculate confidence
    waveform = (trigw-mean(trigw)) / std(trigw); % normalize .wav waveform
    tmp = trigRecord(syncPosition(end)+1:syncPosition(end)+length(waveform)); % find the recorded waveform
    recordWaveform = (tmp-mean(tmp)) / std(tmp); % normalize recorded waveform
    
    
    confidence = abs((recordWaveform*waveform) / (waveform'*waveform)); % calculate correlation between normed .wav and recorded waveform
    if troubleshoot
        plot(waveform); hold on; plot(recordWaveform); hold off;
    end
    
    
    if confidence > thresh
        
        if troubleshoot
            disp(confidence)
            plot(waveform); hold on; plot(recordWaveform); legend('from file','from recording');
            hold off;
        end
        
        % Adjust manually if necessary
        if manualcheck
            
            readytomoveon = 0; firsttimethru = 1;
            while ~readytomoveon
                
                % Plot and select a better alignment location
                figure;
                plot(waveform); hold on; plot(recordWaveform); legend('from file','from recording'); title(confidence);
                hold off;
                if firsttimethru
                    readytomoveon = input('Good to move on? ');
                    firsttimethru = 0;
                end
                
                if ~readytomoveon
                    [newSyncPoint,~] = ginput(1);
                    % Update the syncpoint given this input
                    syncPosition(end) = ceil(syncPosition(end)+newSyncPoint);
                    lastSyncPoint = syncPosition(end); % carry the chosen syncPoint forward
                    lastLength = floor(length(trigw)/20); % step the length forward 20 (to be used if the confidence isn't good enough)
                    
                    % Replot and calculate confidence
                    waveform = (trigw-mean(trigw)) / std(trigw); % normalize .wav waveform
                    tmp = trigRecord(syncPosition(end)+1:syncPosition(end)+length(waveform)); % find the recorded waveform
                    recordWaveform = (tmp-mean(tmp)) / std(tmp); % normalize recorded waveform
                    confidence = abs((recordWaveform*waveform) / (waveform'*waveform)); % calculate correlation between normed .wav and recorded waveform
                    figure;
                    plot(waveform); hold on; plot(recordWaveform); legend('from file','from recording'); title(confidence);
                    hold off;
                    
                    % Check if this was a good enough adjustment
                    readytomoveon = input('Good to move on? ');
                end
            end
            clear readytomoveon;
        else
            if troubleshoot
            disp(confidence);
            end
        end
        
        % Update flag and report result (warn if used Hilbert)
        trigfoundlocs = [trigfoundlocs lastSyncPoint];
        disp(['Confidence ' num2str(confidence)]);
        
        % Plot current progress
        plot(waveform); hold on; plot(recordWaveform); hold off;
        lastLength =length(trigw); % set last length to length of .wav
        wholewaveform(lastSyncPoint:lastSyncPoint+length(waveform)-1) = waveform;
        
    end
    allconfidence = [allconfidence confidence];

    
    
end

singletrigger = zeros(size(trigRecord));
singletrigger(trigfoundlocs) = 1;

figure('color','w');
wholerecord = (trigRecord-mean(trigRecord)) / std(trigRecord); % normalize recorded waveform
plot(wholerecord,'b'); hold on; plot(wholewaveform,'r'); hold on; plot(singletrigger,'k','LineWidth',3);
title('Record: Blue      From File: Red');





end

