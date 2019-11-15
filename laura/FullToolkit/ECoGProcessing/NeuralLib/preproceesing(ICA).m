%% make event for single blocks

Blocks={'B1','B2','B3','B4','B5','B6'};
for i=1:length(Blocks)
    DataPath='./';
    Subject='CUEEG3';
    Block{1}=Blocks{i};
    SoundPath='~/Documents/MATLAB/Nimalab/CUSounds';
    evnt=NeuralFindEvent (DataPath, SoundPath, Subject, Block);
    save(['evnt_',Subject,'_',Blocks{i},'.mat'], 'evnt')
end

%% High pass filtering the data (0.5Hz to higher)
elects=1:62;
cond='htkraw';
new_cond='0.5_higher';
datapath='./';
Blcks={'B1','B2','B3','B4','B5','B6'};
for cnt=1:length(Blcks)
    Blck=Blcks{cnt};
    mkdir([datapath,Blck,'/',new_cond]);
    for cnt1 = elects
        [signal,f] = readhtk([datapath,Blck,'/' cond '/Ch' num2str(cnt1),'.htk']);
        if f~=2400
            display('ERROR');
        end
        
        %%%% WRITE THE FUNCTION HERE %%%%
        signal=signal(:)-mean(signal);
        signalnew=resample(signal,1,24);
        f=100;
        filtered_data=eegfilt(signalnew',f,0.5,0,0,100);
        data_out=filtered_data;
        %%% END Of Function
        
        writehtk([datapath,Blck,'/',new_cond,'/Ch', num2str(cnt1),'.htk'],data_out,f);
        
    end
end


%% Making the out structure for the whole dataset
% datapath='./';
% Subject='CUEEG3';
% cond={'0.5_higher'};
% load(['evnt_',Subject,'.mat']);
% elects=1:62;
% befaft=[0.5 0.5];
% dataf=100;
% specflag='Auditory';
% datatype='EEG';
% %find_artifact();
% out=NeuralGenOut(evnt, datapath, cond, elects, befaft, dataf, specflag, datatype,1);
% save(['./update1/out_',Subject,'_',cond{1},'.mat'], 'out')

% %% Running Kurtosis algorithm on the whole data set
% % run kurtosis three times
% [ALLEEG, EEG, CURRENTSET] =eeglab;
% origin_folder='update1';
% update_folder='preprocessing1';
%
% OutToEEGLab(origin_folder,update_folder,['out_',Subject,'_',cond{1},'.mat']);
%
% %% Save bad Channels
% EEGLabToOut('preprocessing1');
%
% %% find the bad channel locations from update folder
%
% tmp=load('./update1/chanlocs.mat');
% chanloc_origin=tmp.chanlocs;
% tmp=load('./preprocessing1/chanlocs.mat');
% chanloc_pr=tmp.chanlocs;
%
% good_channels=[];
% for i=1:length(chanloc_origin)
%     for j=1:length(chanloc_pr)
%         if strcmp(chanloc_origin(i).labels,chanloc_pr(j).labels)
%             good_channels=[good_channels i];
%             break
%         end
%     end
% end
%%
good_channels=1:60;
%% Make new out structures (for each block separately), with new electrodes
Blocks={'B1','B2','B3','B4','B5','B6'};
datapath='./';
Subject='CUEEG3';
cond={'0.5_higher'};

for i=1:length(Blocks)
    load(['evnt_',Subject,'_',Blocks{i},'.mat']);
    elects=good_channels;
    befaft=[0.5 0.5];
    dataf=100;
    specflag='Auditory';
    datatype='EEG';
    %find_artifact();
    out=NeuralGenOut(evnt, datapath, cond, elects, befaft, dataf, specflag, datatype,1);
    save(['out_',Subject,'_',cond{1},'_',Blocks{i},'.mat'], 'out')
end

%% Run ICA for each block separately
block=6;
[ALLEEG, EEG, CURRENTSET] =eeglab;
cond='0.5_higher';
Blocks={'B1','B2','B3','B4','B5','B6'};
origin_folder='update1_ICA';
update_folder='update2_ICA';
OutToEEGLab(origin_folder,update_folder,['out_',Subject,'_',cond,'_',Blocks{block},'.mat']);


%% Remove ICA componants based on the number you see in previouse step
W_ICA=EEG.icawinv;
remove_comp=[1 5 10 15];
componants=1:60;
componants_keep=setdiff(componants,remove_comp);
cond={'0.5_higher'};
elects=1:60;%good_channels;
Blck=Blocks{block};
cnt1=0;
tmp=[];
for cntelec = elects;
    cnt1=cnt1+1;
    for cnt2=1:length(cond)
        [tmp(:,cnt1),f] = readhtk([datapath,Blck,'/' cond{1} '/Ch' num2str(cntelec)]);
    end
end
signal=tmp;
figure();
topoplot(EEG.icawinv(:,remove_comp(1)),EEG.chanlocs);

ICA=(EEG.icaweights(componants_keep,:)*EEG.icasphere)*signal.';
ICA_purned=EEG.icawinv(:, componants_keep)*ICA;
figure();
plot(signal(:,1));
hold on;
plot(ICA_purned(1,:),'r');

%% Save the new signal in a new directory
f=100;
new_cond='ICA_purned';
datapath='./';
mkdir([datapath,Blck,'/',new_cond]);
for cnt1 = elects
    data_out=ICA_purned(cnt1,:);
    writehtk([datapath,Blck,'/',new_cond,'/Ch', num2str(cnt1)],data_out,f);
end

%% Generate OUT structure
Blocks={'B1','B2','B3','B4','B5','B6'};
datapath='./';
Subject='CUEEG3';
cond={'ICA_purned'};
load(['evnt_',Subject,'.mat']);
elects=1:60;
befaft=[0.5 0.5];
dataf=100;
specflag='Auditory';
datatype='EEG';
%find_artifact();
out=NeuralGenOut(evnt, datapath, cond, elects, befaft, dataf, specflag, datatype,1);
save(['./process1/out_',Subject,'_',cond{1},'.mat'], 'out')
chanlocs=EEG.chanlocs;
save('./process1/chanlocs.mat','chanlocs');

%% find bad channels bad trials
% Running Kurtosis algorithm on the whole data set
% run kurtosis three times
[ALLEEG, EEG, CURRENTSET] =eeglab;
origin_folder='process1';
update_folder='process2';

OutToEEGLab(origin_folder,update_folder,['out_',Subject,'_',cond{1},'.mat']);

%% save reject marks
reject=false(1,size(EEG.data,3));
rej_cond={'rejjp';'rejkurt';'rejmanual';'rejthresh';'rejconst';'rejfreq'};

for i=1:length(rej_cond)
    if ~isempty(EEG.reject.(rej_cond{i}))
        reject=reject|EEG.reject.(rej_cond{i});
    end
end
save(['./' update_folder '/reject'],'reject');

%% save bad chaanels
badchannels=[1 3 4 7 22 29 34 ];
good_channels=setdiff(1:60,badchannels);
save(['./' update_folder '/good_channels'],'good_channels');
%% savechanlocs
chanlocs=EEG.chanlocs;
save(['./' update_folder '/chanlocs'],'chanlocs');

%% finding artifact
artif_name='2';
find_artifact(artif_name);
%%
elects=1:60;
cond='ICA_purned';
new_cond='ICA_purned_normalized';
datapath='./';
Blcks={'B1','B2','B3','B4','B5','B6'};
m=zeros(6,60);
s=zeros(6,60);
for cnt=1:length(Blcks)
    Blck=Blcks{cnt};
    mkdir([datapath,Blck,'/',new_cond]);
    [artifact,fr] = readhtk([datapath,Blck,'/artifact/',artif_name]);
    artifact(artifact==1)=nan;
    for cnt1 = elects
        [signal,f] = readhtk([datapath,Blck,'/' cond '/Ch' num2str(cnt1)]);
        if round(f)~=100
            display('ERROR');
        end
        
        %%%% WRITE THE FUNCTION HERE %%%%
        signalnew=signal(:);
        % signalnew=resample(signal,1,24);
        %         filtered_data=eegfilt(signalnew',100,2,20);
        
        signalnew2=signalnew+artifact(:);
        
        m1=nanmean(signalnew2);
        s1=nanstd(signalnew2);
        
        m(cnt,cnt1)=nanmean(signalnew2);
        s(cnt,cnt1)=nanstd(signalnew2);
        
        %data_out=(signalnew-m1)/s1;
        %%% END Of Function
        
        % writehtk([datapath,Blck,'/',new_cond,'/Ch', num2str(cnt1)],data_out(:)',f);
        
    end
end

% % %%
s_block=mean(s);
for cnt=1:length(Blcks)
    Blck=Blcks{cnt};
    
    
    for cnt1 = elects
        [signal,f] = readhtk([datapath,Blck,'/' cond '/Ch' num2str(cnt1)]);
        if round(f)~=100
            display(['ERROR Frequency sampling is ' num2str(f)]);
        end
        
        %%%% WRITE THE FUNCTION HERE %%%%
        signalnew=((signal-m(cnt,cnt1))/s(cnt,cnt1));
        
        data_out=signalnew;
        %%% END Of Function
        
        writehtk([datapath,Blck,'/',new_cond,'/Ch', num2str(cnt1)],data_out(:)',f);
        
    end
end

%% Generate different bands and their amplitude
elects=1:60;
cond='ICA_purned_normalized';
new_cond_cell={'hilbert_delta','hilbert_theta','hilbert_alpha','hilbert_beta'};
Freq_bands={[1 4],[4 8],[8 15],[15 30]};
for j=1:length(new_cond_cell)
    
    new_cond=new_cond_cell{j};
    freqRange=Freq_bands{j};
    datapath='./';
    Blcks={'B1','B2','B3','B4','B5','B6'};
    for cnt=1:length(Blcks)
        Blck=Blcks{cnt};
        mkdir([datapath,Blck,'/',new_cond]);
        for cnt1 = elects
            [signal,f] = readhtk([datapath,Blck,'/' cond '/Ch' num2str(cnt1)]);
            f=round(f);
            if f~=100
                display('ERROR');
            end
            
            %%%% WRITE THE FUNCTION HERE %%%%
            
            NumOfBands=1;
            doNotch=1;
            dh2 = NeuralHilbert (signal,f,freqRange,NumOfBands,doNotch);
            data_out=dh2(:)';
            %%% END Of Function
            
            writehtk([datapath,Blck,'/',new_cond,'/Ch', num2str(cnt1)],data_out,f);
            
        end
    end
end




% %% find the bad channel locations from update folder
% 
% tmp=load('./update1/chanlocs.mat');
% chanloc_origin=tmp.chanlocs;
% tmp=load('./preprocessing1/chanlocs.mat');
% chanloc_pr=tmp.chanlocs;
% 
% good_channels=[];
% for i=1:length(chanloc_origin)
%     for j=1:length(chanloc_pr)
%         if strcmp(chanloc_origin(i).labels,chanloc_pr(j).labels)
%             good_channels=[good_channels i];
%             break
%         end
%     end
% end


%% making out structures removing bad channels and bad trials and saving them % load dataset saved in process2 first

Blocks={'B1','B2','B3','B4','B5','B6'};
load(['evnt_',Subject,'.mat']);
datapath='./';
Subject='CUEEG3';
conds={'hilbert_theta','hilbert_alpha','hilbert_beta'};%'hilbert_delta'};
load(['./process2/good_channels']);
load(['./process2/reject']);
load(['./process2/reject']);
load(['./process2/mapping_EEGLab.mat'])
cond={};
for cnt=1:length(conds)
    cond{1}=conds{cnt}
dat=EEG.data;
dat(:,:,reject)=[];

% mapping(reject)=[];

elects=good_channels;
befaft=[0.5 0.5];
dataf=100;
specflag='Auditory';
datatype='EEG';
out=NeuralGenOut(evnt, datapath, cond, elects, befaft, dataf, specflag, datatype,1);
out2=out;

data=[];
cnt=1;
for i=1:length(out)
    data(:,1:size(out(i).resp,2),cnt:cnt+size(out(i).resp,3)-1)=out(i).resp;
    cnt=cnt+size(out(i).resp,3);
end

for i=1:length(out)
    out(i).resp=[];
end

mapping2=mapping;


data(:,:,reject)=[];
mapping2(reject)=[];


for i=1:length(out)
    ind=find(mapping2==i);
    if isempty(ind)
        continue;
    end
    if abs(size(out(i).aud,2)-length(find(dat(1,:,ind(1)))))
        abs(size(out(i).aud,2)-length(find(dat(1,:,ind(1)))))
        display('ERROR');
    end
    out(i).resp=data(:,1:size(out(i).aud,2),ind);
end

ind=[];
for i=1:length(out)
    if isempty(out(i).resp)
        display(['block' int2str(i) 'removed']);
        ind=[ind i];
    end
end
out(ind)=[];

save(['out_',Subject,'_',cond{1},'.mat'], 'out')
end

%%
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

%%
 [ALLEEG, EEG, CURRENTSET] =eeglab;
origin_folder='process3';
update_folder='process4';

OutToEEGLab(origin_folder,update_folder,['out_',Subject,'_',cond{1},'.mat']);
