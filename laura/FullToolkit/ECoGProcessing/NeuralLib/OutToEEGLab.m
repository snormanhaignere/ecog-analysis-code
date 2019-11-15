function []=OutToEEGLab(origin_folder,update_folder,out_name)
% Specify origin and update before running this code
% call EEGLab before running this code
% origign_folder= you read the out from
% update_folder= you write the out in
% out_name : name of the out structure
% example: [ALLEEG, EEG, CURRENTSET] =eeglab;
%           OutToEEGLab(origin_folder,update_folder,'out_EEGlab.mat')
% By Bahar 2015

if ~exist('origin_folder') || isempty(update_folder)
    path='./';
else
    path=['./' origin_folder];
    
end

if ~exist('out_name') || isempty(out_name)
   out_name = 'out_EEGlab.mat';
end

load([path '/' out_name])

global EEG
global ALLEEG
global CURRENTSET

[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);

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


display(['Data Sampling: ' int2str(out(1).dataf)]);
display(['Time Points Per Epoch: ' int2str(len)]);
display(['Number of Channels: ' int2str(size(out(1).resp,1))]);

mkdir(['./' update_folder]);
save(['./' update_folder '/Data_EEGLab'],'data');
save(['./' update_folder '/Mapping_EEGLab'],'mapping');
save(['./' update_folder '/Out_EEGLab'],'out');

EEG.setname=update_folder;
EEG.nbchan=channels;
EEG.srate=out(1).dataf;
EEG.trials=epochs;
EEG.pnts=len;
EEG.xmin=0-out(1).befaft(1);
EEG.data=data;

%
%load([path '/chanlocs']);
%EEG.chanlocs=chanlocs;
CURRENTSET=1;
%
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
eeglab redraw;
end
