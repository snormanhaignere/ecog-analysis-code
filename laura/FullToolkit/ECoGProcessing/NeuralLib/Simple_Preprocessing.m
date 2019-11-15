%% HTK 
ConverToHTK(2400);

%% make event for single blocks ( To sync the audio files with recorded audio by analog box) 

Blocks={'B1','B2','B3','B4','B5'};
DataPath='./';
Subject='CUEEG14';
SoundPath='~/Documents/MATLAB/Nimalab/CUSounds';
evnt=NeuralFindEvent (DataPath, SoundPath, Subject, Blocks);
save(['evnt_',Subject,'.mat'], 'evnt')


%% filtering the data 
elects=1:62;
cond='htkraw';
new_cond='2_15hz';
datapath='./';
Blcks={'B1','B2','B3','B4','B5','B6'};
for cnt=1:length(Blcks)
    Blck=Blcks{cnt};
    mkdir([datapath,Blck,'/',new_cond]);
    for cnt1 = elects
        [signal,f] = readhtk([datapath,Blck,'/' cond '/Ch' num2str(cnt1) '.htk']);
        if f~=2400
            error('ERROR');
        end
        
        %          Wn =0.1/500;                   % Normalozed cutoff frequency
        %         [B,A] = butter(4,Wn,'high');
        %         signal=filtfilt(B,A,signal);
        
        %%%% WRITE THE FUNCTION HERE %%%%
        signal=(signal(:)-mean(signal))/std(signal(500:end));
        signalnew=signal;
        signalnew=resample(signal,1,24);
        f=100;
        filtered_data=eegfilt(signalnew',f,2,15,0,100);
        data_out=filtered_data;
        %%% END Of Function
        
        writehtk([datapath,Blck,'/',new_cond,'/Ch', num2str(cnt1) '.htk'],data_out,f);
        
    end
end

%%
% 
datapath='./';
Subject='CUEEG14';
cond={new_cond};%,'hilbert5-10_normalized','hilbert10-20_normalized','hilbert20-35_normalized'};
load(['evnt_',Subject,'.mat']);
elects=1:62;
befaft=[0.5 0.5];
dataf=100;
specflag='Auditory';
datatype='EEG';

 %find_artifact('1');



out=NeuralGenOut(evnt, datapath, cond, elects, befaft, dataf, specflag, datatype);


  for i=1:length(out)
      
      if find(ismember(out(i).name,'K'))
          i
          out(i).name(find(ismember(out(i).name,'K')))='k';
      end
  end
save(['out_',Subject,'_',cond{1}], 'out')  
  
