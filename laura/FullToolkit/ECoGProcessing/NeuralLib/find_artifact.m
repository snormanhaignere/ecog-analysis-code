function []=find_artifact(artifact_name,datapath,condition,Blcks,overlap,step,thrshhold)
%% this function finds the artifact in EEG signal based on a comparison between std of the EEG
% channels and a threshhold, the easy use of this function go to the
% subject directory and run find_artifact;

if ~exist('datapath') || isempty(datapath)
    datapath='./';
end

if ~exist('Blcks') || isempty(Blcks)
    Blcks={'B1','B2','B3','B4','B5','B6'};
end
    


if ~exist('overlap') || isempty(overlap)
    overlap=0.7;
end

if ~exist('stp') || isempty(stp)
    step=0.4;
end

if ~exist('thrshhold') || isempty(thrshhold)
    thrshhold=1.5;
end

if ~exist('condition') || isempty(condition)
    condition='ICA_purned';
end



Blck_name='artifact';

for cnt=1:length(Blcks)
    Blck=Blcks{cnt};
    mkdir([datapath,Blck,'/',Blck_name]);
    [y,f]=readhtk([datapath,Blck,'/',condition,'/Ch1']);
    stp=step*f;
    artifact=zeros(length(y),1);
    base=1;
    
    while(1)
        if((base+stp)>length(y))
            break
        end
        mean_window=mean(y(base:base+stp));
        %artifact(base:base+step)=artifact(base:base+step)|(abs(y(base:base+step)-repmat(mean_window,size(y(base:base+step),1),size(y(base:base+step),2)))>thrshold);
        x=std(y(base:base+stp))>thrshhold*std(y(200:end));
        artifact(base:base+stp)=artifact(base:base+stp)|x;
        base=fix((1-overlap)*stp)+base;
        artifact=artifact(:)';
    end
    writehtk([datapath,Blck,'/artifact/' artifact_name],artifact,2400);
    figure();
    plot(y/20)
    hold on
    plot(artifact*20,'r')
end