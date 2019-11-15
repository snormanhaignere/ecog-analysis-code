function [] = ConvertToHTK_Laura(datapath, datatype, blocks, channels, dataf)
%       Converts ECoG data of various formats to standard HTK format.
%
% Inputs:
%       datapath (optional): subject folder
%           datapath must contain raw data in an 'original' subdirectory
%               eg: datapath / original <- contains raw ECoG
%           default: current directory
%       datatype: specify data format
%           options: 'tdt' 'blackrock' 'xltek' 'natus' 'eeg'
%       blocks (optional): specify which blocks to process
%           format: array of integers OR cell array of block folder names
%               eg: [1:4] OR {'B1' 'B2' 'B3' 'B4'}
%           default: all blocks in datapath / original
%       channels (optional): specify which channels to process
%           format: array of integers
%           default: to all channels in each block
%       dataf (optional): specify desired data fs; best if a multiple of 1000
%           default: 1000Hz
%
% Result:
%       Creates 'processed' directory in patient file.
%       Creates directories for each block.
%       Reads and resamples audio and ECoG data from datapath/original directory.
%       Creates .htk files for each audio and ECoG channel within block file.
%
%       Resulting file structure:
%           Data: datapath / processed / B# / htkraw / Ch#.htk
%           Audio: datapath / processed / B# / audio / Ch#.htk
%
% Hard-Coded Assumptions:
%       1. Audio fs = 24000
%       2. All audio channels will be converted
%       3. Default data channels are all channels in the first block
%       4. Channel and block names can be found by reading files as if they're TDT
%
% Written: 10/2015 by Laura Long
% Last Updated: 11/10/2015 by Laura
% Corrected audfs in "hardcoded assumptions" (was 2400, should be 24000)
% Corrected writehtk block for audio- previously was rewriting a1 each
% channel; added input order to help comments
% Updated: 11/19/2015 by Laura: Added new way of processing blocks (generate blocknames cell array before loop)
% Updated: 11/28/2015 by Laura: Added option to input blocks as cell; started to look at James's edits from last week, need to incorporate them officially and adjust comments
% Updated 12/03/2015 by Laura
% Added James's memory updates: clearing vars, resample as matrix rather than loop
% Added option to loop over TDT channels, rather than the previous hack!
% Added info about inputs to comments
% Updated 9/1/2017 by Laura
% Changed resample mechanisms- now defaults to 1k data 24k audio unless the original sampling rate is lower than desired
% Now all resamples use rat to find factors directly (including TDT, no factor, just calculated directly from desired fs)
% Added option to visualize channels before selecting aud/trig for EDF
% Added option to use WriteChNames_Laura if chnames were not extracted
%
% NOTE: do not modify this function if you are not Laura, for any changes you need to talk to her (nima)

%% Deal with inputs

fprintf('\nConvertToHTK\n');

% Set up original and processed directories

% Check if datapath was entered; if not, default to current directory
if ~exist('datapath','var') || isempty(datapath)
    datapath = pwd; % if not specified, default to the current directory
end

% Define original and destination paths from datapath
origpath = [datapath filesep 'original' filesep];
destpath = [datapath filesep 'processed' filesep];

% Throw error if origpath isn't already there
if ~exist(origpath,'dir')
    error('original folder cannot be found; check datapath');
end

% Create processed file if not already present
if ~exist(destpath,'dir')
    disp('Making processed directory')
    mkdir(destpath);
end

% Set parameter defaults
% Data sampling rate defaults to 1000 and auddefai;tss hard coded to 24000)
if ~exist('dataf','var') || isempty(dataf) % sampling frequency default to 1000
    dataf = 1000;
end
if ~exist('audf','var') || isempty(audf) % sampling frequency default to 1000
    audf = 24000;
end

% Blocks defaults to all blocks found in the original folder; populate blocknames
if ~exist('blocks','var') || isempty(blocks) % blocks default to all
    blocks = getblockinds(origpath); %getblockinds finds vector of block numbers from datapath/original/
end
if iscell(blocks) % If blocks input is matrix, turn it into a cell array
    blocknames = blocks;
else
    for i = 1:length(blocks)
        blocknames{i} = ['B' num2str(blocks(i))];
    end
end

% Channels defaults to all channels found in the first block folder
if ~exist('channels','var') || isempty(channels) % channels default to all
    switch datatype
        case {'tdt1' 'tdt3'}
            tdtchannels = 0;
        case 'gtec'
            channels = 1:63;
        otherwise
            channels = getchannelinds([origpath num2str(blocknames{1})],datatype);
    end
else
    tdtchannels = channels; % hack for TDT for right now- when i figure out how to do the channels individually i'll fix it
end


%% Loop over blocks, reading files and creating new directories for each one
for i = 1:length(blocknames)
    
    % Generate block name, make directories
    disp(['Processing ' blocknames{i} '...']);
    
    % Check whether this block folder already exists; if so, double-check user wants to continue and overwrite
    if exist([destpath, blocknames{i}],'dir')
        processblock = input('This block folder already exists. Overwrite it? ','s');
        if ~strcmpi(processblock,{'y', 'yes', '1'})
            disp('To avoid overwrite, revise block list and rerun');
            break
        end
    end
    
    % Make htkraw directory
    mkdir([destpath, blocknames{i}, filesep, 'htkraw']);
    
    % Check datatype!!! Extract data and set sampling factors accordingly
    switch lower(datatype)
        case 'tdt3'
            % Read all .sev files in the block folder
            disp('...loading .sev files');
            
            
            for j = 1:length(tdtchannels)
                
                % Read the data
                rawdata = SEV2mat([origpath blocknames{i}],'CHANNEL',tdtchannels(j),'VERBOSE',0); % extract data from .sev files in that block
                
                % On the first pass, ID which fields are data vs audio
                if j == 1
                    allfields = fields(rawdata);
                    disp(allfields);
                    datafields = input('Which fields are data? ');
                    datafields = allfields(datafields);
                    audfields = input('Which fields are audio? ');
                    audfields = allfields(audfields);
                end
                
                % Pull out unfiltered ECoG channel and fs for all data fields (each field is an entry in cell array)
                for k = 1:length(datafields)
                    datafield = datafields{k};
                    if isfield(rawdata,datafield) && isfield(rawdata.(datafield),'data');
                        data{k}(j,:) = double(rawdata.(datafield).data); % raw data
                        fs_data(k,j) = double(rawdata.(datafield).fs); % ECoG fs
                    end
                end
                
                % Repeat for audio
                for k = 1:length(audfields)
                    audfield = audfields{k};
                    if isfield(rawdata,audfield) && isfield(rawdata.(audfield),'data')
                        aud{k}(j,:) = double(rawdata.(audfield).data);
                        fs_aud(k,j) = double(rawdata.(audfield).fs);
                    end
                end
                
                clear rawdata;
            end
            
            % Now concatenate the cells together
            datacell = data; data = [];
            for j = 1:length(datacell)
                data = cat(1,data,datacell{j});
            end
            fs_data = unique(fs_data);
            fs_data(fs_data==0) = [];
            if length(fs_data)>1
                disp(fs_data); error('Multiple fs_datas');
            end
            if size(data,1) > length(channels)
                channels = 1:size(data,1);
            end
            
            audcell = aud; aud = [];
            for j = 1:length(audcell)
                aud = cat(1,aud,audcell{j});
            end
            fs_aud = unique(fs_aud);
            if length(fs_aud)>1
                disp(fs_aud); error('Multiple fs_auds');
            end
            
            
        case 'tdt2'
            disp('...loading .tev files');
            
            rawdata = TDTbin2mat([origpath blocknames{i}]);
            
            aud = double(rawdata.streams.Wav5.data);
            fs_aud = double(rawdata.streams.Wav5.fs);
            data = double(rawdata.streams.EEG1.data);
            fs_data = double(rawdata.streams.EEG1.fs);
            if isfield(rawdata.streams,'EEG2');
                try
                    data = cat(1, data, double(rawdata.streams.EEG2.data));
                catch
                    data = cat(1, data, double(rawdata.streams.EEG2.data(1:length(data))));
                end
            end
            
            clear rawdata;
            
        case 'tdt1'
            % Read all .sev files in the block folder
            disp('...loading .sev files');
            
            if tdtchannels == channels % if tdtchannels is the same as input channels, then do it in a loop
                
                for j = 1:length(tdtchannels)
                    rawdata = SEV2mat([origpath blocknames{i}],'CHANNEL',tdtchannels(j),'VERBOSE',0); % extract data from .sev files in that block
                    
                    % Pull out unfiltered ECoG channel and fs
                    data(j,:) = double(rawdata.xWav.data); % raw data
                    fs_data = double(rawdata.xWav.fs); % ECoG fs
                    
                    % Pull out audio channel and fs, if it exists for this channel index
                    if isfield(rawdata,'Aud_') && isfield(rawdata.Aud_,'data')
                        aud(j,:) = double(rawdata.Aud_.data);
                        fs_aud = double(rawdata.Aud_.fs);
                    elseif isfield(rawdata,'Aud1') && isfield(rawdata.Aud1,'data')
                        aud(j,:) = double(rawdata.Aud1.data);
                        fs_aud = double(rawdata.Aud1.fs);
                    elseif isfield(rawdata,'xWv2') && isfield(rawdata.xWv2,'data')
                        aud(j,:) = double(rawdata.xWv2.data);
                        fs_aud = double(rawdata.xWv2.fs);
                    elseif ~isfield(rawdata,'Aud_') && ~isfield(rawdata,'Aud1') && ~isfield(rawdata,'xWv2')
                        error('Cannot find ''Aud_'' or ''Aud1''. Check original file for auditory data.');
                    end
                    
                    clear rawdata;
                end
                
            elseif tdtchannels == 0 % otherwise, if the flag is 0, do all channels at the same time
                
                rawdata = SEV2mat([origpath blocknames{i}],'CHANNEL',tdtchannels,'VERBOSE',0); % extract data from .sev files in that block
                
                % Pull out unfiltered ECoG channel, ECoGfs, audio, audiofs
                if isfield(rawdata,'xWav')
                    data = double(rawdata.xWav.data); % raw data
                    fs_data = double(rawdata.xWav.fs); % ECoG fs
                elseif isfield(rawdata,'ECoG') % or pull filtered, only if xWav doesn't exist (added by Laura 04/26/17 in response to 044_LIJ118, which has only ECoG)
                    data = double(rawdata.ECoG.data); % raw data
                    fs_data = double(rawdata.ECoG.fs); % ECoG fs
                else
                    error('No data field found in rawdata.')
                end
                channels = 1:size(data,1); % if you read all the data, channels should be 1:num channels
                
                if isfield(rawdata,'Aud_')
                    aud = double(rawdata.Aud_.data);
                    fs_aud = double(rawdata.Aud_.fs);
                elseif isfield(rawdata,'Aud1')
                    aud = double(rawdata.Aud1.data);
                    fs_aud = double(rawdata.Aud1.fs);
                elseif isfield(rawdata,'xWv2');
                    aud = double(rawdata.xWv2.data);
                    fs_aud = double(rawdata.xWv2.fs);
                else
                    error('Cannot find ''Aud_'' or ''Aud1''. Check original file for auditory data.');
                end
                
                clear rawdata;
            end
            
            
        case 'edf'
            
            % If I already made the data file, ask whether to load that instead of original files
            if exist([origpath blocknames{i} filesep 'data.mat'],'file')
                loadfromdata = input('Load from data.mat? ');
            else
                loadfromdata = 0;
            end
            
            if loadfromdata
                disp('Loading from data.mat.');
                load([origpath blocknames{i} filesep 'data.mat']); % load data
                
                % Get fs for each var if they don't exist
                if ~exist('fs_data','var'), fs_data = input('What is fs_data? '); end % data
                if exist('aud','var') && ~exist('fs_aud','var'), fs_aud = input('What is fs_aud? '); end % aud
                if exist('trig','var') && ~exist('fs_trig','var'), fs_trig = input('What is fs_trig? '); end % trig
                
            else % Otherwise make from scratch
                
                % Read the datafile
                datafile = dir([origpath blocknames{i} filesep '*.edf']);
                if length(datafile) > 1
                    error('Multiple .edf files. Please align and save as data.mat.');
                end
                datafile = datafile.name;
                [hdr, record] = edfread([origpath blocknames{i} filesep datafile]);
                
                if isfield(hdr,'label')
                    for j = 1:length(hdr.label)
                        disp([num2str(j) hdr.label(j)]);
                    end
                end
                
                % Define data
                datach = input('Which channels are data? ');
                data = record(datach,:);
                try allfs_data = hdr.samples(datach);
                    allfs_data = unique(allfs_data);
                    if length(allfs_data) == 1
                        fs_data = allfs_data;
                    else
                        error('Multiple fs for data channels');
                    end
                catch
                    fs_data = input('What is fs_data? ');
                end
                chnames = hdr.label(datach);
                
                % Option of seeing channels to determine
                visualizetrig = input('Visualize channels to find trigger / audio? ');
                if visualizetrig
                    vischan = input('Which channels do you want to see? ');
                    figure;
                    for cnt = 1:length(vischan)
                        subplot(length(vischan),1,cnt)
                        plot(record(vischan(cnt),:));
                        title(num2str(vischan(cnt)));
                    end
                end
                
                % Define aud and/or trig
                trigoraud = input('Enter 1 for audio, 2 for trigger, or 3 for both: ');
                
                % Get audio if applicable
                if trigoraud == 1 || trigoraud == 3
                    audch = input('Which channels are audio? ');
                    aud = record(audch,:);
                    try allfs_aud = hdr.samples(audch);
                        allfs_aud = unique(allfs_aud);
                        if length(allfs_aud) == 1
                            fs_aud = allfs_aud;
                        else
                            error('Multiple fs for aud channels');
                        end
                    catch
                        fs_aud = input('What is fs_aud? ');
                    end
                end
                
                % Get trigger if applicable
                if trigoraud == 2 || trigoraud == 3
                    trigch = input('Which channels are trigger? ');
                    trig = record(trigch,:);
                    try allfs_trig = hdr.samples(trigch);
                        allfs_trig = unique(allfs_trig);
                        if length(allfs_trig) == 1
                            fs_trig = allfs_trig;
                        else
                            error('Multiple fs for data channels');
                        end
                    catch
                        fs_trig = input('What is fs_trig? ');
                    end
                end
                
                % Save data
                disp('...saving data.mat');
                switch trigoraud
                    case 1
                        save([origpath blocknames{i} filesep 'data.mat'],'data','fs_data','chnames','aud','fs_aud');
                        
                    case 2
                        save([origpath blocknames{i} filesep 'data.mat'],'data','fs_data','chnames','trig','fs_trig');
                        
                    case 3
                        save([origpath blocknames{i} filesep 'data.mat'],'data','fs_data','chnames','aud','fs_aud','trig','fs_trig');
                        
                    otherwise
                        error('Please restart and enter 1, 2, or 3 for trigger vs. audio.')
                end
                
            end
            
            clear record;
            clear hdr;
            
        case 'gtec'
            
            datafile = dir([origpath blocknames{i} filesep '*.hdf5']);
            datafile = datafile.name;
            rawdata = h5read([origpath blocknames{i} filesep datafile], '/RawData/Samples');
            
            data = rawdata(1:63,:);
            data = double(data); %% IS THIS OK?
            fs_data = 4800; %% NEED TO ADJUST THIS; currently inflexible
            
            aud = rawdata(65,:);
            aud = double(aud);
            fs_aud = 4800; %% NEED TO ADJUST THIS; currently inflexible
            
            reference = rawdata(64,:);
            reference = double(reference);
            
            
        case 'blackrock'
            
            % If I already made the data file, load it
            if exist([origpath blocknames{i} filesep 'data.mat'],'file')
                disp('...loading from data.mat');
                load([origpath blocknames{i} filesep 'data.mat']); % load data
                
            else % Otherwise make from scratch
                % Load data
                data = []; aud = [];
                fnames = dir([origpath blocknames{i} filesep '*.ns*']);
                for l = 1:length(fnames)
                    
                    % Read the file
                    disp(['...opening file ' num2str(l) '/' num2str(length(fnames)) ': ' fnames(l).name]);
                    fname = [origpath blocknames{i} filesep fnames(l).name];
                    file = openNSx(fname,'read');
                    
                    % Display channels and ask user which data is aud/data; load accordingly
                    disp(['Comment: ' file.MetaTags.Comment]); disp('Channels: '); disp(cat(2,{file.ElectrodesInfo.ElectrodeID}',{file.ElectrodesInfo.Label}'));
                    audch = input('Which channels are audio? ');
                    datach = input('Which channels are data? ');
                    if ~isempty(audch)
                        fs_aud = file.MetaTags.SamplingFreq;
                        thisaud = [];
                        if iscell(file.Data)
                            for ll = 1:length(file.Data)
                                thisaud = cat(2,thisaud,double(file.Data{ll}(audch,:)));
                            end
                        else
                            thisaud = cat(2,thisaud,double(file.Data(audch,:)));
                        end
                        aud = cat(2,aud,thisaud);
                    end
                    if ~isempty(datach)
                        fs_data = file.MetaTags.SamplingFreq;
                        thisdata = [];
                        if iscell(file.Data)
                            for ll = 1:length(file.Data)
                                thisdata = cat(2,thisdata,double(file.Data{ll}(datach,:)));
                            end
                        else
                            thisdata = cat(2,thisdata,double(file.Data(datach,:)));
                        end
                        data = cat(2,data,thisdata);
                        chnames = {file.ElectrodesInfo(datach).Label};
                    end
                    
                end
                
            end
            
        otherwise
            error(['Datatype "' datatype '" not recognized.']);
            
    end
    
    
    %% Resample and write to HTK
    
    % CORRECTION THAT SHOULD BE USED FOR LIJ126 OR INCORRECT FS IN .SEV ONLY
    %     warning('HACK to fix incorrect fs for LIJ126');
    %     keyboard
    %     load('/Users/LauraLong/Documents/Lab/ECoG Data/057_LIJ126/LIJ126fs.mat')
    %     fs_data = fs;
    
    % Data
    % Resample data unless the orig sampling is less than desired
    if fs_data <= dataf
        dataf = fs_data;
        data_resampled = data;
    else
        disp(['...resampling data to ' num2str(dataf) 'Hz']);
        [upsamp_data,downsamp_data] = rat(dataf/fs_data,1e-10);
        data_resampled = resample(data', upsamp_data, downsamp_data)'; % resample ECoG to dataf using
    end
    clear data;
    
    % Write htk files for each specified data channel
    disp('...writing data to HTK');
    if  ~exist('channels','var') || isempty(channels) % get channels if not populated (applies to blackrock)
        channels = 1:size(data_resampled,1);
    end
    for j = 1:length(channels) % loop to make htk for each channel
        writehtk([destpath, blocknames{i}, filesep, 'htkraw' filesep, 'Ch', int2str(channels(j)), '.htk'], data_resampled(j,:), dataf);
    end
    clear data_resampled;
    
    % Audio files if present
    if exist('aud','var')
        mkdir([destpath, blocknames{i}, filesep, 'analog']);
        
        % Resample audio unless orig sampling rate is lower than desired
        if fs_aud <= audf
            audf = fs_aud;
            aud_resampled = aud;
        else
            disp(['...resampling audio to ' num2str(audf) 'Hz']);
            [upsamp_aud,downsamp_aud] = rat(audf/fs_aud,1e-10);
            aud_resampled = resample(aud', upsamp_aud, downsamp_aud)'; % resample aud to desired freq
        end
        clear aud;
        
        % Write htk files for each audio channel
        disp('...writing audio to HTK');
        for j = 1:size(aud_resampled,1)
            writehtk([destpath, blocknames{i}, filesep, 'analog', filesep, 'a', int2str(j), '.htk'], aud_resampled(j,:), audf); % audio htk
        end
        clear aud_resampled;
        
    end
    
    % Trigger files if present
    if exist('trig','var')
        mkdir([destpath, blocknames{i}, filesep, 'trigger']);
        % Write htk files for each trigger channel
        disp('...writing trigger to HTK');
        for j = 1:size(trig,1)
            writehtk([destpath, blocknames{i}, filesep, 'trigger', filesep, 't', int2str(j), '.htk'], trig(j,:), fs_trig); % trigger htk
        end
    end
    
    % Do reference resample and htk if present
    if exist('reference','var')
        mkdir([destpath, blocknames{i}, filesep, 'reference']);
        disp('...resampling reference');
        reference_resampled = resample(reference', m_data*upsamp_data, downsamp_data)';
        clear reference;
        disp('...writing reference to HTK');
        for j = 1:size(reference_resampled,1)
            writehtk([destpath, blocknames{i}, filesep, 'reference' filesep, 'Ch', int2str(channels(j)), '.htk'], reference_resampled(j,:), dataf);
        end
    end
    
    % Save chnames if present
    if exist('chnames','var') && ~isempty(chnames)
        disp('...saving chnames');
        chnames = strtrim(chnames);
        save([datapath filesep 'original' filesep blocknames{i} filesep 'allchnames_' blocknames{i} '.mat'],'chnames');
    else
        writechnames = input('Write channel names using trode lengths / names? '); % if the chnames weren't loaded before, option to put in the info to the WriteChNames funciton which will concatenate and save to appropriate location
        if writechnames
            trodenames = input('Input trode names as cell array: ');
            trodelens = input('Input trode lengths as array: ');
            WriteChNames_Laura(origpath(1:end-10),blocks,trodenames,trodelens);
        end
    end
    
    disp(['Completed ' blocknames{i}]);
    
    
end

fprintf('Completed ConvertToHTK \n');

end


