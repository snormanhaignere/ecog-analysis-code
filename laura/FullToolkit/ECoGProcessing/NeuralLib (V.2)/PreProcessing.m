%% Bahar 2015
% Generatng EVNT, filtered signal, and out structure

ConvertToHTK(2400)

%%
DataPath='./';
Subject='CUEEG25';
Blocks={'B1','B2','B3','B4'};
SoundPath='/Users/baharkhalighinejad/Desktop/Term3/Inprogress/multispeaker/Task/Sounds';
StimName='StimOrderA.mat';
evnt = NeuralFindEvent (DataPath, SoundPath, Subject, Blocks,[],StimName)
save(['evnt_' Subject '_attended']);

%%

cond1='ms1_howOne_fs1_howTwo';
cond2='fs1_howFour_ms1_howThree';
cond3='ms1_howThree_fs1_howFour' ;
cond4='fs1_howTwo_ms1_howOne';
cond5='ms1_howFour_fs1_howThree';
cond6='fs1_howOne_ms1_howTwo';
cond7='ms1_howTwo_fs1_howOne';
cond8='fs1_howThree_ms1_howFour';






for i=1:length(evnt)
    name=evnt(i).name;
    tmp=strsplit(name,'.wav');
    tmp=strsplit(tmp{1},'_');
    number=tmp{length(tmp)};
    
    if strfind(name,cond1)
        name=['ms1_howOne' number  '_16k.wav'];
    elseif strfind(name,cond2)
        name=['fs1_howFour' number  '_16k.wav'];
    elseif strfind(name,cond3)
        name=['ms1_howThree' number  '_16k.wav'];
    elseif strfind(name,cond4)
        name=['fs1_howTwo' number  '_16k.wav'];
    elseif strfind(name,cond5)
        name=['ms1_howFour' number  '_16k.wav'];
    elseif strfind(name,cond6)
        name=['fs1_howOne' number  '_16k.wav'];
    elseif strfind(name,cond7)
        name=['ms1_howTwo' number  '_16k.wav'];
    elseif strfind(name,cond8)
        name=['fs1_howThree' number  '_16k.wav'];
    end
    evnt(i).name=name;
    evnt(i).StimPath='~/Documents/MATLAB/Nimalab/CUSounds';
end



save(['evnt_' Subject '_attended']);

%%
elects=1:62;
cond='htkraw';
new_cond='2_15hz';
datapath='./';
Blcks={'B1','B2','B3','B4'};
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

datapath='./';
Subject='CUEEG25';
cond={new_cond};%,'hilbert5-10_normalized','hilbert10-20_normalized','hilbert20-35_normalized'};
load(['evnt_',Subject,'_attended.mat']);
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
save(['out_',Subject,'_attended',], 'out')  
