function [phoneme_data,phoneme_feat,senames,dur] = ...
    NeuralPhonemeRep(s, out, phoneme, fname, param,average_single,probabilityint, feats,loc,gender,phnpattern)
% s = ECSpeech;
% out = data (output .mat file from DNN)
% phoneme = phoneme label
% param = [t_before, t_after]
% average_single = 'Average' , 'Single' ('Average' if we have several trial
% of one sentence take the average of all, 'Single' separate trials from
% each other. Default= 'Single' %Bahar_change
% flag = 0, auditory data
% flag = 1, lfp data
% also, combine ?h and ?x and call it ?xh
% gender: if specified, only do m or f subjects.

if ~exist('param','var') || isempty(param)
    t_before = .1;
    t_after  = .6;
else
    t_before = param(1);
    t_after  = param(2);
end

if ~exist('average_single','var') || isempty(average_single)  %Bahar_change
    average_single = 'Single';
end

if ~exist('loc','var') || isempty(loc)
    loc = 'start';
end
if ~exist('feats','var') || isempty(feats)
    feats  = [];
end
if ~exist('gender','var') || isempty(gender)
    gender = [];
end
senames = []; dur = [];
lst = SpeechEvents(s,'phonemes','list');
idx = find(ismember({lst.Note}, phoneme)); % find index of lst with input phoneme
if ~exist('phnpattern','var') || isempty(phnpattern)
    phn = SpeechEvents(s,'phonemes',lst(idx).Note);
    
else
    phn1 = SpeechEvents(s,'phonemes',lst(idx).Note);
    temp=struct2cell(phn1);
    Notes=squeeze(temp(1,1,:));
    temp=struct2cell(phn1);
    Starts=squeeze(temp(3,1,:));
    ind=[];
    for i=1:length(phnpattern)
        
        index=find(ismember(Notes,phnpattern(i).Note));
        for j=1:length(index)
            if phnpattern(i).StartTime==Starts{index(j)}
                ind=[ind index(j)];
            end
        end
        
    end
    phn1(ind)=[];
    phn=phn1;
end

%contains file name, trial number, startTime, and stopTime of files with
%input phoneme
if exist('probabilityint','var') && ~isempty(probabilityint)
    
    tmp=struct2cell(phn);
    tmp=squeeze(tmp(5,1,:));
    prob=zeros(1,length(tmp));
    for i=1:length(tmp)
        if ~isempty(tmp{i})
            prob(i)=tmp{i};
        else
            prob(i)=0;
        end
    end
    
    ind=find(prob>probabilityint(1) & prob<=probabilityint(2));
    phn=phn(ind);
end

srate = out(1).dataf; %sampling rate
t_before = ceil(t_before*srate);
t_after  = ceil(t_after*srate);
phoneme_data = zeros(size(out(1).(fname),1),t_before+t_after,length(phn)); %make_change length(phn) is not correct
%phoneme_data = zeros(size(out(1).(fname)',1),t_before+t_after,length(phn));
% Bahar_change : commented the above line, uncommented the line before that
phoneme_feat = zeros(8,t_before+t_after,length(phn));
onset = ceil(get(s,'PreStimSilence')*srate);
offset = ceil(get(s,'PostStimSilence')*srate);
ind = 1;
names = get(s,'names');
if ~isempty(feats)
    featdata = get(s,feats);
end
for i_phn = 1:length(phn)
    if (isempty(gender)) || (~isempty(gender) && strcmpi(phn(i_phn).Note(1),gender))
        % find the sentence in "out"
        tmp = find(ismember({out.name},phn(i_phn).Note));
        if length(tmp)>1
            warning('not in standard out format');
            idx=tmp(1);
        else
            idx=tmp;
        end
        
        % calculate start index
        start = round(phn(i_phn).StartTime * srate);
        stop  = round(phn(i_phn).StopTime  * srate);
        tmp = round( (start+stop)/2);
        switch loc
            case 'end'
                timepts = [ stop-t_before+1 stop+t_after];
            case 'center',
                timepts = [ tmp-t_before+1 tmp+t_after];
                % capture data window, pad with zeros if problems
            case 'start'
                timepts = [start-t_before+1 start+t_after];
        end
        % if data is missing, record
        if ~isempty(idx)
            if (strcmp(average_single,'Single'))
                trial_numbers=size(out(idx).(fname),3); %Bahar_change
            else
                trial_numbers=1; %Bahar_change
            end
            for idx2=1:trial_numbers %Bahar_change
                data_length =  size(out(idx).(fname),2); %Bahar_change : removed the transpose of out(idx).(fname)'
                if (timepts(1)>onset+5) && (timepts(2)<(data_length-offset+15))
                    
                    
                    % determine which data set to use (flag)
                    % data_to_avg = aud or lfp
                    if (strcmp(average_single,'Average')) %Bahar_change
                        data_to_avg = mean(out(idx).(fname),3); %Bahar_change remove transpose out(idx).(fname)'
                    elseif (strcmp(average_single,'Single')) %Bahar_change
                        data_to_avg =out(idx).(fname)(:,:,idx2);
                    else
                        error(' average_single option is not defined correctly');
                    end
                    
                    
                    data_length =  size(out(idx).(fname),2); %Bahar_change remove transpose out(idx).(fname)'
                    
                    % if the data starting point is before onset, skip:
                    % extract data; if the windows exceed data indices,
                    % pad with zeros
                    if timepts(1)<1 || timepts(2)>data_length
                        continue;
                        %                         if timepts(1)<1 % phoneme too close to start
                        %                             data = data_to_avg(:,1:timepts(2));
                        %                             data = [zeros(size(data,1),abs(timepts(1))+1) data];
                        %                         elseif timepts(2)>data_length % phoneme too close to end
                        %                             data = data_to_avg(:,timepts(1):end);
                        %                             data = [data zeros(size(data,1),timepts(2)-data_length)];
                        %                         end
                    else
                        data = data_to_avg(:,timepts(1):timepts(2));
                        %                         if (trial_numbers==1)
                        %                             artifact=out(idx).artifact(:,:,idx2);
                        % %                             if any(artifact(:,timepts(1):timepts(2)))
                        % %                                 continue;
                        % %                             end
                        %                         end
                        %
                        
                    end
                    phoneme_data(:,:,ind) = data;
                    senames{ind} = phn(i_phn).Note;
                    if ~isempty(feats)
                        senind = find(strcmpi(names,senames{ind}));
                        phoneme_feat(:,:,ind) = featdata{senind}(timepts(1):timepts(2),:)';
                    end
                    dur(ind) = phn(i_phn).StopTime - phn(i_phn).StartTime;
                    ind = ind+1;
                end
            end
        end
    end
end
phoneme_data(:,:,ind:end) = [];
