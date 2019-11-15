%% User Inputs: Set parameters, choose what actions to take

% Enter subject ID.
subject = '025_LIJ110';

% Enter dpath. Should be the patient folder, containing 'original' directory.
% If empty, defaults to current directory.
dpath = ['/Users/Nima/Research/data/eCogProjects2/PhoneRestoration/' subject];
dpath = ['/Users/LauraLong/Documents/Lab/ECoG Data/' subject];

% Enter additional parameters.
% If empty, defaults to all blocks and all channels.
blocks = [55];
channels = [1:64];
window = [];
amplifiers = [];

% Flag to convert to HTK if not already done
converttohtk = 0;

%%  Convert to HTK
if converttohtk == 1
    
    % Set params
    recordingsystem = 'tdt';
    dataf = [];
    
    % Convert
    ConvertToHTK_Laura(dpath,recordingsystem,blocks,channels,dataf);
    
end


%% Generate Report

AnalyzeDataQuality(dpath, subject, blocks, channels, amplifiers, window);