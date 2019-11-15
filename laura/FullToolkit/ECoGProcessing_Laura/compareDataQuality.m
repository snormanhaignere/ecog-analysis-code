%% Script Info:
% Extracts and combines out structures from a subject for different conditions and blocks
% Finds t-values
% Allows user to select which conditions and channels to investigate further
% Calculates STRFs and average unit responses for kept conditions/channels
% Saves tval, strfs, units + subject info + chosen conditions/channels


%% Set Params

subject = '037_NY610';
conditions = {'cleaned_highgamma'}; {'cleaned1_highgamma' 'cleaned_highgamma' 'highgamma'};
blocks = 2:4;

extractdata = 1;
findtvals = 1;
findstrfs = 0;
findavgphn = 1;
findrecon = 1;

%% Load and concatenate all blocks into one out, and all conditions into 'outs'

if extractdata
    
    outs = cell(1,length(conditions));
    for i = 1:length(conditions)
        out = struct([]);
        for j = 1:length(blocks)
            outname = ['/Users/LauraLong/Documents/Lab/ECoG Data/' subject '/processed/B' num2str(blocks(j)) '/out_' subject '_B' num2str(blocks(j)) '_' conditions{i} '.mat'];
            thisout = load(outname);
            out = cat(2,out,thisout.out);
        end
        out = zscore_out(out); % normalize the data
        out = compressresp_out(out); % compress the data
        outs{i} = out; % store
    end
    
end


%% T Values

if findtvals
    
    % Calculate t values for each condition
    tval = cell(1,length(outs));
    for i = 1:length(outs)
        tval{i} = calcSpeechTVal_out(outs{i});
    end
    
    % Extract bad channels from high gamma if necessary
    if strcmp('highgamma',conditions{end}) % check which kind of variables are saved
        clean = load(['/Users/LauraLong/Documents/Lab/ECoG Data/' subject '/processed/B' num2str(blocks(1)) filesep subject '_dataqualityvars.mat']);
        if isfield(clean,'cleanvars')
            tval{end}(clean.cleanvars.channels.bad) = [];
        else
            tval{end}(clean.badchannels) = [];
        end
    end
    
    % Plot all t-values together
    figure;
    for i = 1:length(tval)
        plot(tval{i},'*-'); hold on;
    end
    xlabel('Electrodes'); ylabel('T Values');
    legend(conditions)
    suptitle(subject)
    
    % Plot t-value distribution (can use to determine which outs to use and which threshold to set
    figure;
    for i = 1:length(tval)
        [a,~] = sort(tval{i},'descend');
        plot(a,'.-'); hold on;
    end
    xlabel('Channel'); ylabel('T Value');
    legend(conditions);
    suptitle('Max to Min T Values');
    
end


%% Select channels to keep working with

% Set parameters
whichouts = 1; % which out structures to continue using (by ind)
threshold = 10; % threshold t-value that will be used
numchan = 6; % minimum number of channels to be calculated

% Sort and select channels
[a,b] = sort(tval{whichouts(1)},'descend'); % sort again, this time by the first kept out structure
if sum(a>threshold) > numchan
    bestchan = b(a>threshold);
else
    bestchan = b(1:numchan);
end

% Create outs to keep analyzing
keptouts = outs(whichouts);


%% STRFs

if findstrfs
    
    % Calculate STRFs for each out
    for i = 1:length(keptouts)
        strfs{i} = calcSTRFs_out(keptouts{i},bestchan,[],1,0,1);
    end
    
    % Plot STRF correlations
    figure;
    for j = 1:length(strfs)
        plot(mean(strfs{j}.r_values),'o-'); hold on;
    end
    xlabel('Electrode'); ylabel('STRF Correlation');
    legend(conditions(whichouts));
    suptitle('STRF Correlations by Electrode');
    
    % Plot STRFs side-by-side for each channel
    for i = 1:length(bestchan)
        figure;
        for j = 1:length(strfs)
            subplot(1,length(strfs),j)
            imagesc(squeeze(mean(strfs{j}.strf(:,:,:,i)))); axis xy
            title({conditions{j} , ['STRF Corr ' num2str(mean(strfs{j}.r_values(:,i)))]});
        end
        colormap jet;
        suptitle({['Channel ' num2str(bestchan(i))], ['tval: ' num2str(tval{1}(bestchan(i)))]});
    end
    
    % Plot each condition's STRFs on one figure
    for j = 1:length(strfs)
        figure;
        subinds = numSubplots(length(bestchan));
        for i = 1:length(bestchan)
            subplot(subinds(1),subinds(2),i)
            imagesc(squeeze(mean(strfs{j}.strf(:,:,:,i)))); axis xy
            title({['ch ' num2str(bestchan(i))], ['tval: ' num2str(a(i))] ['corr ' num2str(mean(strfs{j}.r_values(:,i)))]})
        end
        colormap jet
        suptitle(conditions{j});
    end
    
end

%% Average Phoneme Responses

if findavgphn
    
    % Calculate average phoneme responses; will also plot average spectrograms and average responses by phoneme
    units = cell(1,length(whichouts));
    for i = 1:length(whichouts)
        [units{i},s] = calcAvgUnitResp(outs{whichouts(i)},'Phonemes',{'seg13'},[],0,1,0,1,1);
    end
    
    % Plot average responses by channel
    for j = 1:length(units)
        figure;
        subinds = numSubplots(length(bestchan));
        for i = 1:length(bestchan)
            subplot(subinds(1),subinds(2),i)
            im(squeeze(units{j}.mean_resp(bestchan(i),:,:))');
            set(gca, 'YTick', [1:1:length(units{j}.units)], 'YTickLabel', units{j}.units)
            title({['ch ' num2str(bestchan(i))], ['tval: ' num2str(tval{1}(bestchan(i)))]})
        end
        suptitle(['Average Phoneme Responses by Channel: ' conditions{j}]);
        colormap jet
    end
    
    % Plot average compressed responses by channel
    for j = 1:length(units)
        figure;
        subinds = numSubplots(length(bestchan));
        for i = 1:length(bestchan)
            subplot(subinds(1),subinds(2),i)
            im(squeeze(units{j}.mean_resp_compressed(bestchan(i),:,:))');
            set(gca, 'YTick', [1:1:length(units{j}.units)], 'YTickLabel', units{j}.units)
            title({['ch ' num2str(bestchan(i))], ['tval: ' num2str(tval{1}(bestchan(i)))]})
        end
        suptitle(['Average Phoneme Responses by Channel after Compression: ' conditions{j}]);
        colormap jet
    end
    
end

%% Reconstruction

if findrecon
    
    % Run reconstruction
    for i = 1:length(keptouts)
        [recon{i},reconparams{i}] = calcReconstruction_out(keptouts{i},bestchan,[],[],[],[],[],[],[],[],1);
    end
    
    % Plot reconstruction vs original
    whichepochs = 1:5;
    for i = 1:length(recon)
        
        for j = 1:length(whichepochs);
            figure;
            subplot(2,1,1);
            imagesc(recon{i}.stim_tests{whichepochs(j)}); axis xy
            title('Original');
            subplot(2,1,2);
            imagesc(recon{i}.rstims{whichepochs(j)}); axis xy
            xlabel('Time (s)'); ylabel('Frequency Band');
            title(['Reconstruction     Correlation: ' num2str(recon{i}.recons(whichepochs(j)))]);
            xlabel('Time (s)'); ylabel('Frequency Band');
            colormap(jet);
        end
         
    end
    
    
end



%% Save Data if Desired

saveinput = input('Save data? ');
% saveinput = 1;
if saveinput == 1
    
    sizeunit = whos('units');
    if sizeunit.bytes > 1e+9
        strffilename = ['/Users/LauraLong/Documents/MATLAB/NimaLab/Projects/SegmentationDataAnalysis' filesep subject '_tval_strf.mat'];
        unitfilename = ['/Users/LauraLong/Documents/MATLAB/NimaLab/Projects/SegmentationDataAnalysis' filesep subject '_unit.mat'];
        save(strffilename,'subject','conditions','blocks','tval','strfs','whichouts','bestchan');
        disp(['Saved ' strffilename]);
        save(unitfilename,'subject','conditions','blocks','whichouts','bestchan','units','-v7.3');
        disp(['Saved ' unitfilename]);
    else
        filename = ['/Users/LauraLong/Documents/MATLAB/NimaLab/Projects/SegmentationDataAnalysis' filesep subject '_tval_strf_unit.mat'];
        save(filename,'subject','conditions','blocks','tval','strfs','whichouts','bestchan','units');
        disp(['Saved ' filename]);
    end
    
    
end

