function[]= convert_hdf5_to_htk(dataf,varargin)
% Convering EEG, ECoG or MEG data to standard HTK format
% Usage: can be used in two ways:
%       1: ConvertToHtk ( fs, original path, destination path)
%       original path: contains the "original" folder of raw data, assuming the standard naming
%       scheme (B1, B2, ...)
%       2: Go to "CUEEGxx" foler which contains the "original" folder and
%       run  ConvertToHTK([]) or ConvertToHTK(fs)
% dataf is the sampling frequency of data, default is 2400 Hz
% Raw data is saved in y(channel*time) variable, y(2:63,:) is the EEG data,
% y(64:65,:) are the references, y(66,:) is the analog channel
% 
% modified from ConvertToHTK.m by Tasha ??
% May 2015

% Setting the Original and destination defaults
org='./';
dest='./';

% Setting the default Frequency Sampling value
if ~exist('dataf') || isempty(dataf)
    dataf=2400;
end

if nargin==2
    error('Enter original and destination folder')
elseif nargin>2
    org=[varargin{1},'/'];
    dest=[varargin{2},'/'];
end

%reading all files in "original" folder
files=dir([org,'original']);

%extracting the names of the files
filenames=cell(length(files),1);
for i=1:length(files)
    filenames{i}=files(i).name;
end

%extracting data files
cnt=1;
DataFiles=[];
for i=1:length(filenames)
    if length(filenames{i}) > 4 && strcmp(filenames{i}(end-3:end), 'hdf5')
        DataFiles{cnt}=filenames{i};
        cnt=cnt+1;
    end
end


nametmp=' ';
cnt=1;
signal=[];


for i=1:length(DataFiles)
    % checking if the current file is a segment of the previous file:
    if ~isempty(strfind(nametmp,DataFiles{i}(1:end-5))) 
        cnt=cnt+1;
        
        %check if we have a loaded data from the previous step
    elseif ~isempty(signal)
        % the previous step was the last segment of a block, save the block
        % in htk format
        
        mkdir([dest,'B',int2str(part),'/htkraw']);
        mkdir([dest,'B',int2str(part),'/reference']);
        mkdir([dest,'B',int2str(part),'/analog']);
        mkdir([dest,'B',int2str(part),'/Stimulus']);
        writehtk([dest,'B',int2str(part),'/analog/a1.htk'],signal(65,:),dataf);
        for j=1:62
            writehtk([dest,'B',int2str(part),'/htkraw/Ch',int2str(j),'.htk'],signal(j,:),dataf);
        end
        
        writehtk([dest,'B',int2str(part),'/reference/r1.htk'],signal(63,:),dataf);
        writehtk([dest,'B',int2str(part),'/reference/r2.htk'],signal(64,:),dataf);
        display(['Block ' int2str(part) ' generated ...']);
        cnt=1;
        signal=[];
    end
    %find the number of current block
    splitted=strsplit(DataFiles{i}, {'B','_'});
    part=str2double(splitted(2));
    
    %read the data of current block and save it to signal
    filename=[org,'original/',DataFiles{i}];
    signal=h5read(filename,'/RawData/Samples');
    nametmp=DataFiles{i}(1:end-5);
    
    if i==length(DataFiles)
        mkdir([dest,'B',int2str(part),'/htkraw']);
        mkdir([dest,'B',int2str(part),'/reference']);
        mkdir([dest,'B',int2str(part),'/analog']);
        mkdir([dest,'B',int2str(part),'/Stimulus']);
        writehtk([dest,'B',int2str(part),'/analog/a1.htk'],signal(65,:),dataf);
        for j=1:62
            writehtk([dest,'B',int2str(part),'/htkraw/Ch',int2str(j),'.htk'],signal(j,:),dataf);
        end
        
        writehtk([dest,'B',int2str(part),'/reference/r1.htk'],signal(63,:),dataf);
        writehtk([dest,'B',int2str(part),'/reference/r2.htk'],signal(64,:),dataf);
        display(['Block ' int2str(part) ' generated ...']);
        cnt=1;
        signal=[];
    end
    
end

% Finding the stimorders and putting them in blocks
for i=1:length(filenames)
    if ~isempty(strfind(filenames{i},'stimorder'))
        splitted=strsplit(filenames{i}, {'B','_'});
        part=str2double(splitted(2));
        copyfile([org,'original/',filenames{i}],[dest,'B',int2str(part),'/Stimulus/StimOrder.mat']);
    end
end

end