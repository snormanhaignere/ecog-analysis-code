function ConvertXltekToHTK(neuralmatfile,NumOfBlocks,destpath)
% convert Xltek ECoG continous recording to standard blocks
% Note: Xltek uses sync channel to find event times
%
% Inputs: 
%       neuralmatfile:      neuralmat structure with preprocessed neural data
%                           structure should have the following elements:
%       neuralmat.patientID         internal ID, i.e., 'CU0xx'
%       neuralmat.resp              neural response
%                                   dimensions: (N electrodes x time)
%                                   1st electrode = sync channel
%       neuralmat.blocksep          indices separating blocks
%                                   should have dimension NumOfBlocks + 1
%       neuralmat.names             ordered names of sound files
%       neuralmat.channels          channel labels
%       neuralmat.units             channel units, i.e., uV
%       neuralmat.Fs                sampling frequency used in DAQ
%                                   for Xltek, should be 500 Hz
%
%       NumOfBlocks:        number of blocks used (1 for each separate
%                           task)
%       destpath:           destination where .htk files should be saved
%
% Tasha: April 13, 2015

%
%addpath('/Users/tashanagamine/Research/ECoG/stimulus')
%addpath('/Users/tashanagamine/Dropbox/Nimalab/GeneralCode/NeuralLib')
%addpath('/Users/tashanagamine/Dropbox/Nimalab/GeneralCode/PreprocessingGUI/subFunctions')

% function inputs
%neuralmatfile = 'NeuralMat_CU033.mat';
%NumOfBlocks = 4; 
%destpath = pwd;

% load the neural data
load(neuralmatfile);
dataf = neuralmat.Fs;

% number of channels
numchannels = length(neuralmat.channels);

% neural response, size = (N electrodes) x (time)
neuraldataAll = neuralmat.resp; 

% load indices separating blocks
b = neuralmat.blocksep; 

%zeros(1,length(double(neuralmat.(['C' num2str(neuralmat.Header.ChannelID(1))]))) );

for cnt1 = 1:numchannels % channel number
    
    %neuraldata = double(neuralmat.(['C' num2str(neuralmat.Header.ChannelID(cnt1))]));
    neuraldata = neuraldataAll(cnt1,:); 
    display('loaded channel');
    
    for cnt2 = 1:NumOfBlocks % block number
        
        % cut the sound and neural data for each block
        %aud_tmp  = audiodata(tmp(cnt2)*soundf:tmp(cnt2+1)*soundf);
        data_tmp = neuraldata(b(cnt2):b(cnt2+1));
        
        % where to write this?
        datapath = [destpath filesep 'B' num2str(cnt2) filesep 'htkraw' filesep];
        if ~exist(datapath,'dir'), mkdir(datapath);end
        
        %stimpath = [destpath filesep 'B' num2str(cnt2) filesep 'analog' filesep];
        %if ~exist(stimpath,'dir'), mkdir(stimpath);end
        % now write the sound of this block to the "aoppropriate directory
        %writehtk([stimpath filesep 'a1.htk'],aud_tmp,soundf);
        
        htkname = [datapath 'Ch' num2str(cnt1) '.htk'];
        writehtk(htkname,data_tmp,dataf);
        
    end
end

