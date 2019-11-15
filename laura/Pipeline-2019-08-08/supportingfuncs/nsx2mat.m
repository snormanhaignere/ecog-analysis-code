function nsx2mat(varargin)

% Converts NSx files to MAT files and saves as *.nsxmat within the
% specified directory. The function takes up to two inputs. One of the
% inputs is the full filename to the NSx file to be converted. Local disk
% space triple the size of the file is required for this function to run.
% Works for filespec 2.1, 2.2 and 2.3. If you want to downsample, add the
% downsample frequency scalar (i.e. 2000) in units S/sec as an
% additional input. 
%
% Example: nsx2mat('I:\sample.ns5',10000);
% Example: nsx2mat('2000','I:\sample.ns4');
% Example: nsx2mat('I:\sample.ns4');
% Example: nsx2mat(1000);
% Example: nsx2mat;
%
% Version Date: 20131119
%   - update: Elliot Smith (20131119)
% Author: Tyler Davis

FileName = '';
dsFs = []; %downsample frequency
nsx2mat_version = 20131119;

switch nargin    
    case 2
        if isnumeric(varargin{1})
            FileName = varargin{2};
            dsFs = varargin{1};
        else
            FileName = varargin{1};
            dsFs = varargin{2};
        end            
    case 1
        if isnumeric(varargin{1})
            dsFs = varargin{1};
        else
            FileName = varargin{1};
        end       
end

if isempty(FileName)
    [filename,path] = uigetfile('I:\Data\*.ns*','Choose nsx file...');
    FileName = fullfile(path,filename);
end

NSxFullName = FileName;
NSxMatFullName = [NSxFullName,'mat'];
TMPFullName = [NSxFullName(1:end-3),'tmp'];

switch NSxFullName(end-2:end)
    case 'ns5'
        Fs = 30000;
    case 'ns4'
        Fs = 10000;
    case 'ns3'
        Fs = 2000;
    case 'ns2'
        Fs = 1000;
    otherwise
        disp('Choose an NSx file')
        return
end

% Version and sample frequency checking
try
    if exist(NSxMatFullName,'file')==2
        load(NSxMatFullName,'-mat','Header')
        if isfield(Header,'nsx2mat_version')
            if Header.nsx2mat_version~=nsx2mat_version
                delete(NSxMatFullName);
            elseif ~isempty(dsFs)
                if Header.Fs~=dsFs
                    delete(NSxMatFullName);
                end
            else
                if Header.Fs~=Fs
                    delete(NSxMatFullName);
                end
            end
        else
            delete(NSxMatFullName);
        end
    end
    clear('Header')
catch ME
    disp('Cannot delete old nsxmat file. Run Matlab as administrator')
    return
end

% Converting to mat file if all prior conditions are satisfied
if exist(NSxMatFullName,'file')==0
    
    Header.nsx2mat_version = nsx2mat_version;
    
    % Creating filter object
    if ~isempty(dsFs)
        dsFs = Fs/round(Fs/dsFs);
        if Fs > dsFs
            N = 10; %order
            Fpass = dsFs/4; %passband frequency
            Apass = 1; %passband ripple (dB)
            Astop = 80; %stopband attenuation (dB)
            h = fdesign.lowpass('N,Fp,Ap,Ast',N,Fpass,Apass,Astop,Fs);
            Hd = design(h,'ellip');
        else
            disp('Cannot upsample! Check the sampling rate.')
            return
        end
    end
    
    % Checking filespec
    FID = fopen(NSxFullName, 'r', 'l');
    fseek(FID, 8, 'bof');
    filespec = fread(FID, [1,2],   '*uchar');
    fseek(FID, 0, 'bof');
    
    % Reading NSx file    
    if all(filespec==[2,2]) || all(filespec==[2,3])
        
        % Reading filespec 2.2 or 2.3
        Header.FileID       = fread(FID, [1,8],   '*char');
        Header.FileSpec     = fread(FID, [1,2],   '*uchar');
        Header.HeaderBytes  = fread(FID, [1,1],   '*uint32');
        Header.Fs           = fread(FID, [1,16],  '*char');
        Header.Comment      = fread(FID, [1,256], '*char');
        Header.Period       = fread(FID, [1,1],   '*uint32');
        Header.Resolution   = fread(FID, [1,1],   '*uint32');
        Header.TimeOrigin   = fread(FID, [1,8],   '*uint16');
        Header.ChannelCount = fread(FID, [1,1],   'uint32=>double');
        
        for k = 1:Header.ChannelCount
            Header.Type(k,:)           = fread(FID, [1,2],  '*char');
            Header.ChannelID(k,:)      = fread(FID, [1,1],  '*uint16');
            Header.ChannelLabel(k,:)   = fread(FID, [1,16], '*char');
            Header.PhysConnector(k,:)  = fread(FID, [1,1],  '*uint8');
            Header.ConnectorPin(k,:)   = fread(FID, [1,1],  '*uint8');
            Header.MinDigVal(k,:)      = fread(FID, [1,1],  '*int16');
            Header.MaxDigVal(k,:)      = fread(FID, [1,1],  '*int16');
            Header.MinAnlgVal(k,:)     = fread(FID, [1,1],  '*int16');
            Header.MaxAnlgVal(k,:)     = fread(FID, [1,1],  '*int16');
            Header.Units(k,:)          = fread(FID, [1,16], '*char');
            Header.HighFreqCorner(k,:) = fread(FID, [1,1],  '*uint32');
            Header.HighFreqOrder(k,:)  = fread(FID, [1,1],  '*uint32');
            Header.HighFiltType(k,:)   = fread(FID, [1,1],  '*uint16');
            Header.LowFreqCorner(k,:)  = fread(FID, [1,1],  '*uint32');
            Header.LowFreqOrder(k,:)   = fread(FID, [1,1],  '*uint32');
            Header.LowFiltType(k,:)    = fread(FID, [1,1],  '*uint16');
        end
        
        Header.DataHeader           = fread(FID, [1,1], '*uint8');
        Header.DataTimestamp        = fread(FID, [1,1], '*uint32');
        Header.ChannelLengthSamples = fread(FID, [1,1], '*uint32');
        
    else
        
        % Reading filespec 2.1
        Header.FileID       = fread(FID, [1,8],  '*char');
        Header.Label        = fread(FID, [1,16], '*char');
        Header.Period       = fread(FID, [1,1],  '*uint32');
        Header.ChannelCount = fread(FID, [1,1],  'uint32=>double');
        Header.ChannelID    = fread(FID, [Header.ChannelCount,1], '*uint32');
        
    end
    
    BegOfData = ftell(FID);
    fseek(FID, 0, 'eof');
    EndOfData = ftell(FID);
    fseek(FID, BegOfData, 'bof');
    
    DataLengthBytes = EndOfData - BegOfData;
    
    % Determining system memory to maximize data segments
    [~,sV] = memory;
    SystemMemory = sV.PhysicalMemory.Available;
    
%     SystemMemory = regexp(evalc('feature memstats'),'\d*(?= MB)','match');
%     SystemMemory = str2double(SystemMemory{2})*1e6; % Units bytes
    
    % Calculating maximum data segment to load into memory
    SegmentCount = ceil(DataLengthBytes/(0.75*SystemMemory));
    SegmentSamples = round(double(Header.ChannelLengthSamples)/SegmentCount);
    SegmentDivisor = floor(double(Header.ChannelLengthSamples)/SegmentSamples);
    SegmentRemainder = rem(double(Header.ChannelLengthSamples),SegmentSamples);
    if SegmentRemainder==0
        SegmentMatrix = repmat(SegmentSamples,SegmentDivisor,1);
    else
        SegmentMatrix = [repmat(SegmentSamples,SegmentDivisor,1);SegmentRemainder];
    end
    
    % Parsing and saving to *.tmp file
    save(TMPFullName,'Header','-v7.3')
    for k = 1:size(SegmentMatrix,1)
        clc, fprintf('NSX2MAT Parsing: %0.0f%% complete\nLoading %0.1f GB of data\n',((k-1)*Header.ChannelCount*100)/(size(SegmentMatrix,1)*Header.ChannelCount),SegmentMatrix(k)*2/1e7)
        tempData = fread(FID, [Header.ChannelCount,SegmentMatrix(k)], '*int16');
        for m = 1:Header.ChannelCount
            clc, fprintf('NSX2MAT Parsing: %0.0f%% complete\n',(((k-1)*Header.ChannelCount+m)*100)/(size(SegmentMatrix,1)*Header.ChannelCount))
            tempSubChan = ['C',num2str(Header.ChannelID(m)),'_',num2str(k)];
            eval([tempSubChan,' = tempData(m,:);']);
            save(TMPFullName,tempSubChan,'-append','-v7.3')
            clear(tempSubChan)
        end
        clear('tempData')
    end
    
    % Reconstructing from *.tmp file
    if ~isempty(dsFs)
        Header.ChannelLengthSamples = (Header.ChannelLengthSamples-rem(double(Header.ChannelLengthSamples)-1,round(Fs/dsFs))-1)/round(Fs/dsFs)+1;
        Header.Fs = dsFs;
        Header.Period = round(double(Header.Resolution)/dsFs);
    else
        Header.Fs = Fs;
    end
    save(NSxMatFullName,'Header','-v7.3')
    for m = 1:Header.ChannelCount
        clc, fprintf('NSX2MAT Reconstructing: %0.0f%% complete\n',m*100/double(Header.ChannelCount))
        tempChan = ['C',num2str(Header.ChannelID(m))];
        eval([tempChan,' = [];'])
        for k = 1:size(SegmentMatrix,1)
            tempSubChan = ['C',num2str(Header.ChannelID(m)),'_',num2str(k)];
            load(TMPFullName,tempSubChan,'-mat')
            eval([tempChan,' = [',tempChan,',',tempSubChan,'];'])
            clear(tempSubChan)
        end
        %%%%%%%%%%%%%%%%%%%%%%%% Downsample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(dsFs)
            eval([tempChan,' = filtfilt(Hd.sosMatrix,Hd.ScaleValues,double(',tempChan,'));'])
            eval([tempChan,' = ',tempChan,'(1:',num2str(round(Fs/dsFs)),':end);'])
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save(NSxMatFullName,tempChan,'-append','-v7.3')
        clear(tempChan)
    end
    
    delete(TMPFullName)
    fclose(FID);
    
else
    clc, fprintf('nsxmat file already exists!\n')
end



