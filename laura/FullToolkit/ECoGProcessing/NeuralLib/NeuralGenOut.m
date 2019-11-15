function out = NeuralGenOut(evnt, datapath,cond, elects, befaft, dataf, specflag, datatype,artifact)
% evnt is the event strcutre from NeuralFindEvents
% cond is the data condition, e.g. raw, afterICA, ...
% elects is the electrode numbers to be used
% befaft: how much of the data before and after the stimulus to be included
% dataf: out data sampling frequency
% specflag: what type of spectrogram, options: Auditory
% artifact: artifacts that you want to add and are included in artifact
% folder

if ~exist('datapath') || isempty(befaft)
    datapath='./';
end
if ~exist('artifact')
    artifact=[];
end

if ~exist('befaft') || isempty(befaft)
    befaft=[0.5,0.5];
end

if ~exist('elects') || isempty(elects)
    elects = 1:128;
end
if ~exist('datatype','var') || ~isempty(datatype)
    datatype = 'ECoG';
end

if ~exist('dataf','var') || isempty(dataf)
    no_resample=1;
else
    no_resample=0;
end

names=cell(length(evnt),1);
for i=1:length(evnt)
    names{i}=evnt(i).name;
end
% [x,unique_mat,unique_index]=unique(names);
% out=struct('name',cell(1,length(unique_mat)));
out=struct('name',cell(1,length(names)));

P_Blck='  ';
loadload;close;

for cnt=1:length(evnt)
    
%     i=unique_index(cnt);
    i=cnt;
    if isempty(out(i).name)
        str=evnt(cnt).name;
        disp(['Processing sound ',num2str(cnt),': ',str]);
        if strcmp(str(end-2:end),'wav')
            soundchar=strsplit(str, {'/','.wav'});
            out(i).name=soundchar{end-1};
            [audio,audiof] = audioread([evnt(cnt).StimPath evnt(cnt).name]);
        else
            out(i).name = str;
            [audio,audiof] = audioread([evnt(cnt).StimPath '/' evnt(cnt).name,'.wav']);
        end
        
        audio = audio(:,1); 
       
        audio=[zeros(befaft(1)*audiof,1);audio;zeros(befaft(2)*audiof,1)];
        out(i).sound=audio;
        out(i).soundf=audiof;
        out(i).dataf=dataf;
        out(i).duration=evnt(cnt).stopTime-evnt(cnt).startTime+sum(befaft);
        out(i).befaft=befaft;
        out(i).type= datatype;
        out(i).resp=[];
        out(i).artifact=[];
        out(i).channelnames=evnt(i).channelnames;
        
        switch specflag
            case 'Auditory'
                tmpaud = wav2aud(audio, [1000/dataf 1000/dataf -2 log2(audiof/16000)] )';
                
                try
                out(i).aud=tmpaud(:,1:round(out(i).duration*dataf));
                catch
                    warning('Size error on line 79')
                    out(i).aud=tmpaud; % Change James: round(out(i).duration*dataf) was longer than tmpaud
                end
            case 'mfcc'
                out(i).mfcc=melfcc_wrapper(audio,audiof);
        end
        if abs(size(tmpaud,2)-round(out(1).duration*dataf))>1,
            warning('Size mismatch');
        end
        
        try
            out(i).aud=tmpaud(:,1:round(out(i).duration*dataf));
        catch
            warning('Size error on line 79')
            out(i).aud=tmpaud; % Change James: round(out(i).duration*dataf) was longer than tmpaud
        end
        out(i).trial=cnt;
    else
        out(i).trial=[out(i).trial cnt];
    end
    
    Blck=evnt(cnt).block;
    
    if Blck==P_Blck;
        
    else
        tmp = [];
        for cnt1 = 1:length(elects) % edit by Tasha; gets rid of indexing issue
            [tmp(:,cnt1),f] = readhtk([datapath Blck '/' char(cond) '/Ch' num2str(elects(cnt1)) '.htk']);
        end
        if no_resample
            record=tmp;
            dataf=f;
        else
            record=resample(tmp,dataf,round(f));
        end
        %%find_artifact
        tmp=[];
        for cnt1 = artifact           
            [tmp(:,cnt1),f_art] = readhtk([datapath Blck '/artifact/' num2str(cnt1)]);
        end
        if ~isempty(artifact)
            artif_wave=downsample(tmp,round(f_art/dataf));
        end
    end
    
    P_Blck=Blck;
    
    eegstart=round((evnt(cnt).startTime-befaft(1))*dataf);
    newdata=record(eegstart+1:eegstart+round(out(i).duration*dataf),:)';
    out(i).resp=cat(3,out(i).resp,newdata);
    if ~isempty(artifact)
        newartifact=artif_wave(eegstart+1:eegstart+round(out(i).duration*dataf),:)';
        out(i).artifact=cat(3,out(i).artifact,newartifact);
    end
    
    
end
