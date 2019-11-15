

datapath = '/Users/LauraLong/Documents/Lab/ECoG Data/034_LIJ113';
blocks = [216 217]
audchannels = 1:4

fromhtk = 1;
fromsev = 0;

if fromhtk
    
    for i = blocks
        
        figure;
        
        for j = audchannels
            StimPath = [datapath filesep 'Processed' filesep 'B' num2str(i) filesep 'analog' filesep 'a' num2str(j) '.htk']; % define path to audio htk file
            [audioRecord, audioRecordFreq] = readhtk(StimPath); % reads the audio waveform and fs from the htk file
            audioRecord = audioRecord(:)'; % ensure the audio is a row vector
            
            subplot(length(audchannels),1,j);
            plot(audioRecord);
            title(['aud channel ' num2str(j)]);
        end
        
        suptitle(['Audio for Block ' num2str(i)]);
        
    end
    
end

clear aud;

if fromsev
    
    figure;
    
    for j = audchannels
        rawdata = SEV2mat([datapath filesep 'original' filesep 'B' num2str(i)],'CHANNEL',j,'VERBOSE',0); % extract data from .sev files in that block
        
        % Pull out audio channel and fs, if it exists for this channel index
        if isfield(rawdata,'Aud_') && isfield(rawdata.Aud_,'data')
            aud = double(rawdata.Aud_.data);
            fs_aud = double(rawdata.Aud_.fs);
        elseif isfield(rawdata,'Aud1') && isfield(rawdata.Aud1,'data')
            aud = double(rawdata.Aud1.data);
            fs_aud = double(rawdata.Aud1.fs);
        elseif isfield(rawdata,'xWv2') && isfield(rawdata.xWv2,'data')
            aud = double(rawdata.xWv2.data);
            fs_aud = double(rawdata.xWv2.fs);
        elseif ~isfield(rawdata,'Aud_') && ~isfield(rawdata,'Aud1') && ~isfield(rawdata,'xWv2')
            error('Cannot find ''Aud_'', ''Aud1'', or ''xWv2''. Check original file for auditory data.');
        end
        
        clear rawdata;
        
        subplot(length(audchannels),1,j);
        plot(aud);
        title(['aud channel ' num2str(j)]);
    end
    
    suptitle(['Audio for Block ' num2str(i) ': from .sev']);
    
end