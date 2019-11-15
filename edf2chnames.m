function chnames = edf2chnames(exp, subjid, r, varargin)

% Prints and saves the channel names for a given subject's EDF file
% 
% 2019-11-11: Commented, Sam NH

global root_directory

I.print = true;
I = parse_optInputs_keyvalue(varargin, I);

% directory for this project
project_directory = [root_directory '/' exp];

% just read in labels
edf_file = [project_directory '/data/ECoG-EDF/' subjid '/r' num2str(r) '.edf'];
if ~exist(edf_file, 'file')
    error('EDF file does not exist: %s', edf_file);
end
hdr = edfread(edf_file);
chnames = hdr.label;

% save to MAT file
MAT_file = [project_directory '/data/ECoG-chnames/' subjid '.mat'];
save(mkpdir(MAT_file), 'chnames');

% save to text file
TXT_file = [project_directory '/data/ECoG-chnames/' subjid '.txt'];
fid = fopen(mkpdir(TXT_file), 'w');
for i = 1:length(chnames)
    if I.print
        fprintf(fid, '%d: %s\n', i, chnames{i});
        fprintf('%d: %s\n', i, chnames{i});
    end
end
fclose(fid);



