function [rawdata, rawref, rawampref, hgdata, hgref, hgampref] = AnalyzeDataQuality(datapath, subject, blocks, channels, amplifiers, window)
%% Deal with inputs

fprintf('\nAnalyzeDataQuality\n');

% Set parameter defaults

% Check if datapath was entered; if not, default to current directory
if ~exist('datapath','var') || isempty(datapath)
    datapath = pwd; % if not specified, default to the current directory
end

% Blocks defaults to all blocks found in the original folder; populate blocknames
if ~exist('blocks','var') || isempty(blocks) % blocks default to all
    blocks = getblockinds([datapath filesep 'Processed']); %getblockinds finds vector of block numbers from datapath/original/
end
if iscell(blocks) % If blocks input is matrix, turn it into a cell array
    blocknames = blocks;
else
    for i = 1:length(blocks)
        blocknames{i} = ['B' num2str(blocks(i))];
    end
end

% Set ampref to 1 if amplifier input is populated; will reference each amplifier separately
if ~exist('amplifiers','var') || isempty(amplifiers)
    ampref = 0;
else
    ampref = 1;
end



%% Begin looping over blocks, get directories and extract high gamma if necessary

for i = 1:length(blocks)
    
    % Define datapaths
    rawdatapath = [datapath filesep 'Processed' filesep blocknames{i} filesep 'htkraw'];
    hgdatapath = [datapath filesep 'Processed' filesep blocknames{i} filesep 'highgamma'];
    
    % Default channels to all in raw folder
    if ~exist('channels') || isempty(channels) % channels default to all
        channels = getchannelinds(rawdatapath,'htk');
    end
    
    % Extract high gamma if not already done
    if ~exist(hgdatapath,'dir') || isempty(hgdatapath)
        HTKtoHighGamma(datapath,'htkraw',blocks,channels,100);
    end
    
    
    %% Load data
    
    rawdata = [];
    hgdata = [];
    for j = 1:length(channels)
        disp(['Ch' num2str(j)]);
        if ~exist('window','var') || isempty(window)
            [rawdata(j,:),rawfs] = readhtk([rawdatapath filesep 'Ch' num2str(channels(j)) '.htk']);
            [hgdata(j,:),hgfs] = readhtk([hgdatapath filesep 'Ch' num2str(channels(j)) '.htk']);
        else    
            [rawdata(j,:),rawfs] = readhtk([rawdatapath filesep 'Ch' num2str(channels(j)) '.htk'],[window(1) window(2)]*1000);
            [hgdata(j,:),hgfs] = readhtk([hgdatapath filesep 'Ch' num2str(channels(j)) '.htk'],[window(1) window(2)]*1000);
        end
    end
    
    
    % Common Reference with first PC
    rawstd = mapstd(rawdata);
    [u,s] = svd(rawstd*rawstd');
    rawref = rawstd - u(:,1:4)*(rawstd'*u(:,1:4))';
    
    hgstd = mapstd(hgdata);
    [u,s] = svd(hgstd*hgstd');
    hgref = hgstd - u(:,1:4)*(hgstd'*u(:,1:4))';
    
    % Option to Reference each amp separately
    if ampref
        rawampref = zeros(size(rawstd));
        hgampref = zeros(size(hgstd));
        for j = 1:length(amplifiers)
            
            rawamp = rawstd(amplifiers{j},:);
            [u,s] = svd(rawamp*rawamp');
            rawampref(amplifiers{j},:) = rawamp - u(:,1) * (rawamp'*u(:,1))';
            
            hgamp = hgstd(amplifiers{j},:);
            [u,s] = svd(hgamp*hgamp');
            hgampref(amplifiers{j},:) = hgamp - u(:,1) * (hgamp'*u(:,1))';
        end
        
        % Plot Raw and HG with and without References
        fig = figure;
        subplot(6,1,1);
        im(tanh(mapstd(rawdata)));
        xlabel('Time');
        ylabel('Channels');
        title('Raw Data');
        
        subplot(6,1,2);
        im(tanh(rawref));
        xlabel('Time');
        ylabel('Channels');
        title('Raw Data: Referenced with PC1');
        
        subplot(6,1,3);
        im(tanh(rawampref));
        xlabel('Time');
        ylabel('Channels');
        title('Raw Data: Referenced with PC1 by Amplifier');
        
        subplot(6,1,4);
        im(tanh(hgstd));
        xlabel('Time');
        ylabel('Channels');
        title('High Gamma');
        
        subplot(6,1,5);
        im(tanh(hgref(:,:)));
        xlabel('Time');
        ylabel('Channels');
        title('High Gamma: Referenced with PC1');
        
        subplot(6,1,6);
        im(tanh(hgampref));
        xlabel('Time');
        ylabel('Channels');
        title('High Gamma: Referenced with PC1 by Amplifier');
        
        suptitle([subject ': ' blocknames{i}]);
        % fix the underscore problem, think about adding other info to the title
        
    else
        rawampref = [];
        hgampref = [];
        
        % Plot Raw and HG with and without References
        fig = figure;
        subplot(4,1,1);
        % data_rr = resample(rawdata',1,10)';
        im(tanh(mapstd(rawdata)));
        xlabel('Time');
        ylabel('Channels');
        title('Raw Data');
        
        subplot(4,1,2);
        im(tanh(hgref(:,:)));
        xlabel('Time');
        ylabel('Channels');
        title('Raw Data: Referenced with PC1');
        
        subplot(4,1,3);
        im(tanh(hgstd));
        xlabel('Time');
        ylabel('Channels');
        title('High Gamma');
        
        subplot(4,1,4);
        im(tanh(hgref(:,:)));
        xlabel('Time');
        ylabel('Channels');
        title('High Gamma: Referenced with PC1');
        
        suptitle([subject ': ' blocknames{i}]);
        % fix the underscore problem, think about adding other info to the title
        
        
    end

   % Save figure to PDF   
    saveevnt = input('Save figure to PDF? ','s');
    if ~strcmpi(saveevnt,{'y', 'yes', '1'})
    else
        print([datapath filesep subject '_' blocknames{i}],'-dpdf');
        % see if there's a way to make it better fill the PDF screen
        disp(['Saved to ' datapath filesep subject '_' blocknames{i} '.pdf']);
    end
    
    
end

