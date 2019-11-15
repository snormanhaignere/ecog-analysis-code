function[] = EEGLabToOut(update_folder)
global EEG
global ALLEEG

load(['./' update_folder '/Out_EEGLab.mat']) 
load(['./' update_folder '/mapping_EEGLab.mat']) 

if length(mapping)~=size(EEG.data,3);
    error('the mapping in the update folder does not match the EEG file');
end

reject=false(1,size(EEG.data,3));
rej_cond={'rejjp';'rejkurt';'rejmanual';'rejthresh';'rejconst';'rejfreq'};

for i=1:length(rej_cond)
    if ~isempty(EEG.reject.(rej_cond{i}))
        reject=reject|EEG.reject.(rej_cond{i});
        EEG.reject.(rej_cond{i})=[];
    end
end

for i=1:length(out)
    out(i).resp=[];
end



EEG.data(:,:,reject)=[];
mapping(reject)=[];


for i=1:length(out)
    ind=find(mapping==i);
    if isempty(ind)
        continue;
    end
    if abs(size(out(i).aud,2)-length(find(EEG.data(1,:,ind(1)))))
        abs(size(out(i).aud,2)-length(find(EEG.data(1,:,ind(1)))))
        display('ERROR');
    end
    out(i).resp=EEG.data(:,1:size(out(i).aud,2),ind);
end

ind=[];
for i=1:length(out)
    if isempty(out(i).resp)
        display(['block' int2str(i) 'removed']);
        ind=[ind i];
    end
end
out(ind)=[];
chanlocs=EEG.chanlocs;

len=0;
epochs=0;
for i=1:length(out)
    len=max(size(out(i).resp,2),len);
    epochs=epochs+size(out(i).resp,3);
end
channels=size(out(i).resp,1);
data=zeros(channels,len,epochs);
mapping=zeros(epochs,1);
cnt=1;
for i=1:length(out)
    data(:,1:size(out(i).resp,2),cnt:cnt+size(out(i).resp,3)-1)=out(i).resp;
    mapping(cnt:cnt+size(out(i).resp,3)-1)=i;
    cnt=cnt+size(out(i).resp,3);
end

EEG.data=data;
EEG.trials=epochs;
EEG.pnts=len;

ALLEEG=EEG;

[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
eeglab redraw;
save(['./' update_folder '/Out_EEGLab'],'out');
save(['./' update_folder '/mapping_EEGLab'],'mapping');
save(['./' update_folder '/chanlocs'],'chanlocs');
end