%
function []=splitdata(file,Block_name,StimOrder,soundpath,recordedfilename,fECOG,datachannel,triggername,fsrecord)
%file='data1.edf';
%%%Block1
% first_soundFile='ms1_howTwo1_16k.wav';
% recordedfilename='audio1.m4a';
% Block_name='B1';

%%%Block2
% first_soundFile='ms1_howThree1_16k.wav';
% recordedfilename='audio1.m4a';
% Block_name='B2';

%%%Block3
% first_soundFile='ms1_howFour1_16k.wav';
% recordedfilename='audio1.m4a';
% Block_name='B3';
% 
% soundpath='/Users/baharkhalighinejad/Desktop/summer/InProgress/ECogTaskstory/Task/Male_Female_Speakers/Sounds';

first_soundFile=StimOrder{1};
first_soundFile=[soundpath '/' first_soundFile];
if fsrecord==44100
fs=4410;  %% experiment sampling frequency
%fECOG=500;
else
    fs=fsrecord;
end

[y,header]=readedf(file,fECOG);
channelnames = cellstr(header.channelname);

triggerchannel=find(strcmp(channelnames,triggername));

% datachannel=find(~strcmp(channelnames,'Event'));
%%%Aljo 
% datachannel=[2:17,19:34,35:63,65:105,107:109];

trigs = y(triggerchannel,:);
 if strcmp(triggername,'TRIG1')
    trigs(find(trigs>0))=0;
    trigs(find(abs(trigs)>10000))=-120;
    y(triggerchannel,:)=trigs;
 end

trigs(find(abs(trigs)>0.15*std(abs(trigs))*mean(abs(trigs))))=0.15*std(abs(trigs))*mean(abs(trigs));
trigs=abs(trigs)>(0.1*std(abs(trigs))*mean(abs(trigs)));

if isempty(datachannel)
    datachannel=1:size(y,1);
    tmp=std(y,[],2);
    datachannel([triggerchannel;find(tmp==0)])=[];
end

%
figure
plot(trigs)
[T1 ~]=ginput(2);

data=y(:,round(T1(1)):round(T1(2)));

data=data(:,1:end-rem(length(data),fECOG));
%data=data.';
%save('data.mat','data');

plot(data(triggerchannel,:));

%% read the audio File

[audiotmp,fsrecord] = audioread(recordedfilename);
audiotmp=audiotmp(1:end-rem(length(audiotmp),fsrecord));
audio=resample(audiotmp,fs,fsrecord);
figure
plot(1:20:length(audio),audio(1:20:end))
[T1 ~]=ginput(2);
audioFile=audio(T1(1):T1(2));
audioFile=audioFile(1:end-rem(length(audioFile),fs));
clear fsrecord audiotmp T1 audio

%% read the first audiofile

audioRecord=audioFile;

aduioRecordFreq=fs;


syncPosition = [];
audioRecord=audioRecord(:)';
stop = size(audioRecord,2);
lastSyncPoint = 0;
lastLength = 0;
goToNext = 1; % condition flag on confidence, if confidence is big enough then go to the next sentence
soundIndex = 1;
goToNext=1;
tmplast=0;
[sound1,fs_audio] = audioread(first_soundFile);
sound1=resample(sound1,fs,fs_audio);
w2=sound1(:,1);

while 1

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
    
    tmp = audioRecord(syncPosition(end)+1:syncPosition(end)+length(waveform) );
    %         tmp = audioRecord(syncPosition(end)+1:min(syncPosition(end)+length(waveform),length(audioRecord)) );
    recordWaveform = (tmp-mean(tmp)) / std(tmp);
    len_w=round(length(recordWaveform));
    confidence = abs((recordWaveform(1:len_w)*waveform(1:len_w)) / (waveform(1:len_w)'*waveform(1:len_w)));
    %         confidence2 = abs((recordWaveform(end-1000:end)*waveform(end-1000:end)) / (waveform(end-1000:end)'*waveform(end-1000:end)));
    %         confidence=(confidence+confidence2)/2;
    if confidence > 0.3
        break
    end    
end

figure('color','w');
plot(waveform); hold on; plot(recordWaveform,'r');
firstfile_start=(syncPosition(end)+1)/fs;



%% read the last audiofile

% 
% audioRecord=audioFile;
% aduioRecordFreq=fs;
% 
% 
% syncPosition = [];
% audioRecord=audioRecord(:)';
% stop = size(audioRecord,2);
% lastSyncPoint = length(audioFile)-fs*600;
% lastLength = 0;
% goToNext = 1; % condition flag on confidence, if confidence is big enough then go to the next sentence
% soundIndex = 1;
% goToNext=1;
% tmplast=0;
% [sound1,fs_audio] = audioread(last_soundFile);
% sound1=resample(sound1,fs,fs_audio);
% w2=sound1(:,1);
% 
% while 1
% 
%     % Cross-correlation
%     windowSize = ceil(length(w2)/aduioRecordFreq)+0.5;
%     audioPart = audioRecord(lastSyncPoint+lastLength+1:min(lastSyncPoint+lastLength+windowSize*aduioRecordFreq, stop)); % sampling rate is 10K, find the xcorr in the next 20s window
%     c = xcorr(w2, audioPart);
%     syncPoint = abs(length(audioPart)-(find(abs(c)==max(abs(c)))));
%     syncPosition = [syncPosition; syncPoint+lastSyncPoint+lastLength];
%     
%     lastSyncPoint = syncPosition(end)-2;
%     while  tmplast>=lastSyncPoint
%         lastSyncPoint=lastSyncPoint+1;
%     end
%     
%     tmplast=lastSyncPoint;
%     lastLength = 1;%ceil(length(w2)/600); %% 60 is step size and it should change based on the length of the waveform
%     
%     % Calculate confidence
%     waveform = (w2-mean(w2)) / std(w2);
%     
%     tmp = audioRecord(syncPosition(end)+1:syncPosition(end)+length(waveform) );
%     %         tmp = audioRecord(syncPosition(end)+1:min(syncPosition(end)+length(waveform),length(audioRecord)) );
%     recordWaveform = (tmp-mean(tmp)) / std(tmp);
%     len_w=round(length(recordWaveform));
%     confidence = abs((recordWaveform(1:len_w)*waveform(1:len_w)) / (waveform(1:len_w)'*waveform(1:len_w)));
%     %         confidence2 = abs((recordWaveform(end-1000:end)*waveform(end-1000:end)) / (waveform(end-1000:end)'*waveform(end-1000:end)));
%     %         confidence=(confidence+confidence2)/2;
% 
%     if confidence > 0.3
%         break
%     end    
% end
% 
% figure('color','w');
% plot(waveform); hold on; plot(recordWaveform,'r');
% 
% lastfile_start=(syncPosition(end)+1)/fs;

%%
figure
trigs=data(triggerchannel,:);
trigs(find(abs(trigs)>0.15*std(abs(trigs))*mean(abs(trigs))))=0.15*std(abs(trigs))*mean(abs(trigs));
trigs=abs(trigs)>(0.1*std(abs(trigs))*mean(abs(trigs)));
plot(trigs(1:fECOG*400));
[X,Y]=ginput(1)
plot(X(1)-fECOG:X(1)+fECOG,trigs(X(1)-fECOG:X(1)+fECOG))
[X,Y]=ginput(1)


syncP=X(1);
syncS=syncP/fECOG;
lengthS=min((length(data)/fECOG)-syncS,(length(audioRecord)/fs)-firstfile_start);

figure('color','w');

audio=audioRecord((firstfile_start-3)*fs:lengthS*fs+(firstfile_start-3)*fs);
data=data(:,(syncS-3)*fECOG:lengthS*fECOG+(syncS-3)*fECOG);
trigs=trigs(:,(syncS-3)*fECOG:lengthS*fECOG+(syncS-3)*fECOG);

plot(trigs);
hold on


tmpaudio=resample(audio,floor(fs/fECOG)*fECOG,fs);
tmpaudio=downsample(tmpaudio,floor(fs/fECOG));

plot(tmpaudio,'r');

%lastfile_start2=lastfile_start*fECOG-tmpsync+3*fECOG+1;
%line([lastfile_start2 lastfile_start2],[-1 1]);
line([3*fECOG 3*fECOG],[-1 1]);

trigs=data(triggerchannel,:);
trigs(find(abs(trigs)>0.15*std(abs(trigs))*mean(abs(trigs))))=0.15*std(abs(trigs))*mean(abs(trigs));
trigger=abs(trigs)>(0.1*std(abs(trigs))*mean(abs(trigs)));
trigger(1:3*fECOG-50)=0;
triggerind=find(trigger);
removind=find(diff(triggerind)<fECOG);
triggerind(removind+1)=[];
triggertime=triggerind/fECOG;


data=data(datachannel,:);
channelnames=channelnames(datachannel);


save([Block_name '_data_ECoG.mat'],'data');
save([Block_name '_trigger_ECoG.mat'],'triggertime');
save([Block_name '_audio_ECoG.mat'],'audio');
save([Block_name '_channelnames_ECoG.mat'],'channelnames');

end
