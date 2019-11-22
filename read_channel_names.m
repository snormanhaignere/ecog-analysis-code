function chnames = read_channel_names(subjid)

global root_directory;
root_directory = my_root_directory;

speech_TCI_subjid = {'CU106', '086_NY723', 'CU105', 'CU103', '093_CUBF44', '085_CU102', 'CU108', 'NY751', 'CU107'};

if strcmp(subjid, '057_LIJ126')
    chnames = cell(1, 352);
    for i = 1:352
        chnames{i} = ['ch' num2str(i)];
    end
elseif strcmp(subjid, '078_CU100')
    chnames = cell(1, 96);
    for i = 1:96
        chnames{i} = ['ch' num2str(i)];
    end
elseif strcmp(subjid, '080_NY706')
    E = read_elec_coords_pial(subjid, 'both');
    chnames = E.chnames;s
% elseif strcmp(subjid, '093_CUBF44')
%     X = load([root_directory '/speech-TCI/data/ECoG-Raw/' subjid '/allchnames_B1.mat'], 'chnames');
%     chnames = X.chnames;
% elseif strcmp(subjid, 'CU107')
%     X = load([root_directory '/speech-TCI/data/ECoG-Raw/' subjid '/allchnames_B1.mat'], 'chnames');
%     chnames = X.chnames;
% elseif strcmp(subjid, 'CU108')
%     X = load([root_directory '/speech-TCI/data/ECoG-Raw/' subjid '/allchnames_B1.mat'], 'chnames');
%     chnames = X.chnames;
% elseif strcmp(subjid, 'NY751')
%     X = load([root_directory '/speech-TCI/data/ECoG-Raw/' subjid '/allchnames_B1.mat'], 'chnames');
%     chnames = X.chnames;
elseif any(ismember(speech_TCI_subjid, subjid))
    load([root_directory '/speech-TCI/data/ECoG-chnames/' subjid '.mat'], 'chnames');
else    
    X = load([root_directory '/ecog-quilting/data/ECoG-Raw/' subjid '/allchnames_B1.mat'], 'chnames');
    chnames = X.chnames;
end

for i = 1:length(chnames)
    x = cast(chnames{i}, 'uint8');
    x = x(x~=0);
    chnames{i} = cast(x, 'char');
end

channel_names_directory = [root_directory '/ecog-quilting/figures/channel-names'];
if ~exist(channel_names_directory, 'dir'); mkdir(channel_names_directory); end
channel_file = [channel_names_directory '/' subjid '.txt'];
fid = fopen(channel_file, 'w');
for i = 1:length(chnames)
    fprintf(fid, '%4d: %40s\n', i, chnames{i});
end
fclose(fid);

%% Modify to match the coordinate files provided by Hugh

switch subjid
    case '077_NY701'
        
        new_chnames = cell(size(chnames));
        for i = 1:length(chnames)
            
            % number
            num_as_string = regexp(chnames{i}, '(\d)*', 'match');
            num_index = regexp(chnames{i}, '(\d)*', 'start');
            assert(length(num_as_string)==length(num_as_string));
            format_num = strrep(sprintf('%2d',str2double(num_as_string{1})), ' ', '0');
            
            % name
            name = chnames{i}(1:num_index-1);
            if strcmp(name, 'GRID')
                name = 'G';
            end
            new_chnames{i} = [name format_num];
            
        end
        
        chnames = new_chnames;
        
    case '062_NY668'
        
        chnames = strrep(chnames, 'EEG', '');
        chnames = strrep(chnames, 'REF', '');
        chnames = strrep(chnames, '_', '');        
        
    case '065_NY676'
        
        chnames = strrep(chnames, 'EEG', '');
        chnames = strrep(chnames, 'REF', '');
        chnames = strrep(chnames, '_', ''); 
        
        
    case {'068_NY679', '070_NY686', '071_NY682', '074_NY688', '079_NY704', '080_NY712', '081_NY708', '086_NY723'}
        
        new_chnames = cell(size(chnames));
        for i = 1:length(chnames)
            
            % number
            num_as_string = regexp(chnames{i}, '(\d)*', 'match');
            num_index = regexp(chnames{i}, '(\d)*', 'start');
            if isempty(num_as_string)
                new_chnames{i} = chnames{i};
            else
                assert(length(num_as_string)==1);
                format_num = strrep(sprintf('%2d',str2double(num_as_string{1})), ' ', '0');
                
                % name
                name = chnames{i}(1:num_index-1);
                new_chnames{i} = [name format_num];
            end
            
        end
        
        chnames = new_chnames;
                           
end

if strcmp(subjid, '086_NY723')
    chnames = strrep(chnames, 'RDS', 'RDSI');
end

if strcmp(subjid, '080_NY712')
    chnames = strrep(chnames, 'DAM', 'DAMT');
    chnames = strrep(chnames, 'DPM', 'DPMT');
end
    