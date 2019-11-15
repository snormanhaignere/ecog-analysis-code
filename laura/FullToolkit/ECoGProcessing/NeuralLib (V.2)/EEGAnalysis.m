%% Global Params

DataPath = './';
Subject = 'CUEEG2';

%% Evnt

Blocks = {'B1' 'B2' 'B3' 'B4'};
SoundPath = '~/Documents/Lab/EEG Data/Task_Sounds/CUSounds'; 

evnt = NeuralFindEvent(DataPath, SoundPath, Subject, Blocks);
save(['evnt_' Subject '.mat'],'evnt')

%% Out

cond = {'htkraw'};
load(['evnt_' Subject '.mat']);
elects = 1:62;
befaft = [0.5 0.5];
dataf = 100;
specflag = 'Auditory';
datatype = 'EEG';

out = NeuralGenOut(evnt,DataPath,cond,elects,befaft,dataf,specflag,datatype);
save(['out_' Subject '_' cond{1} '.mat'], 'out')





