function out = NeuralGenOut_ECoG(evnt, PatientPath, experiment, elects, befaft, specflag, datatype)
% evnt is the event strcutre from NeuralFindEvents
% cond is the data condition, e.g. raw, afterICA, ...
% elects is the electrode numbers to be used
% befaft: how much of the data before and after the stimulus to be included
% dataf: out data sampling frequency
% specflag (optional): what type of spectrogram, options: Auditory mfcc
% htktype (optional): htkraw / htkfilt - filtered data used for visualization purposes on the TDT.
% resample_flag (optional): use htk data that were resampled to 1000Hz
% artifact: artifacts that you want to add and are included in artifact
% folder

if ~exist('befaft','var') || isempty(befaft)
    befaft=[0.5,0.5];
end
if ~exist('specflag','var') || isempty(specflag)
    specflag='Auditory';
end
if ~exist('elects','var') || isempty(elects)
    elects = 1:128;
end


names=cell(length(evnt),1);
for i=1:length(evnt)
    names{i}=evnt(i).name;
end

fs_ECoG = evnt.fs_ECoG;
fs_aud = evnt.fs_aud;

% [x,unique_mat,unique_index]=unique(names);
% out=struct('name',cell(1,length(unique_mat)));
out=struct('name',cell(1,length(names)));

loadload;close;

P_Blck = ''; % Variable to determine if data needs to be read from HTK

for cnt=1:length(evnt)
    
    i=cnt;
    if isempty(out(i).name)
        str=evnt(cnt).name;
        
        % Read ECoG data from htk if needed
        Blck=evnt(cnt).block; % Get current block
        
        if ~strcmp(Blck,P_Blck); % Determine if data needs to be read
            clear data
            disp('Reading data from HTK')
            idx = 1;
            for elec = elects
                try
                    data(idx,:) = readhtk([PatientPath '/processed/' evnt(i).block '/' datatype '/htkraw/' evnt(i).channelnames{elec} '.htk']);
                catch
                    data(idx,:) = readhtk([PatientPath '/processed/' evnt(i).block '/' datatype '/htkraw/Ch' int2str(elec) '.htk']);
                end
                idx = idx + 1;
            end
            
            P_Blck=Blck; % Reasign P_Blck so that data isn't read until block change            
        end
        
        
        % Segment data
        out(i).duration=evnt(i).stopTime-evnt(i).startTime + sum(befaft);
        resp = data(:,evnt(i).sync - round(befaft(1)*fs_ECoG) : evnt(i).sync - round(befaft(1)*fs_ECoG) + round(out(i).duration*fs_ECoG) -1);

        % Assign data to  out structure
        out(i).trial=cnt;
        out(i).name = evnt(i).name;
        out(i).resp = resp;
        out(i).dataf = fs_ECoG;
        out(i).soundf = fs_aud;
               
        % Get Auditory Spectrograms -------------------------------------
        disp(['Processing sound ',num2str(cnt),': ',str]);
        switch experiment
            case {'Multi_Speaker','Multi_Speaker_2'}
                Stim_Attend = evnt(i).Stim_Attend;
                Stim_Unattend = evnt(i).Stim_Unattend;
                
                Stim_Attend=[zeros(round(befaft(1)*fs_aud),1); Stim_Attend; zeros(round(befaft(2)*fs_aud),1)];
                Stim_Unattend=[zeros(round(befaft(1)*fs_aud),1); Stim_Unattend; zeros(round(befaft(2)*fs_aud),1)];
                
                switch specflag
                    case 'Auditory'
                        tmp_attend = wav2aud(Stim_Attend, [1000/fs_ECoG 1000/fs_ECoG -2 log2(fs_aud/16000)] )';
                        tmp_unattend = wav2aud(Stim_Unattend, [1000/fs_ECoG 1000/fs_ECoG -2 log2(fs_aud/16000)] )';
                        
                        try
                            out(i).aud_attend=tmp_attend(:,1:round(out(i).duration*fs_ECoG));
                            out(i).aud_unattend=tmp_unattend(:,1:round(out(i).duration*fs_ECoG));
                        catch
                            warning('Size error')
                            out(i).aud_attend=tmp_attend; % Change James: round(out(i).duration*dataf) was longer than tmpaud
                            out(i).aud_unattend=tmp_unattend; % Change James: round(out(i).duration*dataf) was longer than tmpaud
                        end
                    case 'mfcc'
                        out(i).mfcc=melfcc_wrapper(audio,fs_aud);
                end
                if abs(size(tmp_attend,2)-round(out(i).duration*fs_ECoG))>1,
                    warning('Size mismatch');
                end
                
                out(i).stim_attend=Stim_Attend;
                out(i).stim_unattend=Stim_Unattend;
                
            otherwise
                stim = evnt(i).stim;
                
                if size(stim,1) > size(stim,2)
                    stim = stim';
                end
                
                % Pad with zeros
                stim=[zeros(1,round(befaft(1)*fs_aud)) stim zeros(1,round(befaft(2)*fs_aud))];
                
                switch specflag
                    case 'Auditory'
                        [tmp] = wav2aud(stim, [1000/fs_ECoG 1000/fs_ECoG -2 log2(fs_aud/16000)] )';
                        
                        try
                            out(i).aud=tmp(:,1:round(out(i).duration*fs_ECoG));
                            % This will usually happen
                            if abs(size(tmp,2)-round(out(i).duration*fs_ECoG))>1,
                                warning('Slight size mismatch - probably fine.');
                            end
                        catch
                            warning('Large size mismatch - double check')
                            out(i).aud=tmp; 
                            out(i).resp = out(i).resp(:,1:size(tmp,2));
                        end
                    case 'mfcc'
                        out(i).mfcc=melfcc_wrapper(audio,fs_aud);
                end
                
                out(i).stim=stim;
        end
        
        % Assign additional data to out structure
        out(i).channelnames = evnt(i).channelnames(elects);
        out(i).befaft=befaft;
        out(i).type = 'ECoG';
        out(i).artifact=[];
        
        clear resp
    else
        out(i).trial=[out(i).trial cnt];
    end
end
