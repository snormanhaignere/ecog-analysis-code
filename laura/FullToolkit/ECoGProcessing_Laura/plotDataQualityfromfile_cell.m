%% Script Info: plots the most important figures from compareDataQuality from the saved param file

%% Choose subject
subject = '039_NY622';
loaddata = 1;
ploteachstrf = 0;
plotavgeachphn = 0;
plotavgspec = 0;
whichphnrespbychan = 'both'; % can be compressed, norm, or both

%% Load data
if loaddata
    filename = ['/Users/LauraLong/Documents/MATLAB/NimaLab/Projects/SegmentationDataAnalysis' filesep subject '_tval_strf_unit.mat'];
    if exist(filename,'file')
        load(filename)
    else
        strffilename = ['/Users/LauraLong/Documents/MATLAB/NimaLab/Projects/SegmentationDataAnalysis' filesep subject '_tval_strf.mat'];
        unitfilename = ['/Users/LauraLong/Documents/MATLAB/NimaLab/Projects/SegmentationDataAnalysis' filesep subject '_unit.mat'];
        load(strffilename)
        load(unitfilename)
    end
end

%% Plot t-values

% Plot original t values
figure;
for i = 1:length(tval)
    plot(tval{i},'*-'); hold on;
end
xlabel('Electrodes'); ylabel('T Values');
legend(conditions)
suptitle([subject ': T Values by Electrode'])

% Plot t-value distribution (can use to determine which outs to use and which threshold to set
figure;
for i = 1:length(tval)
    [a,~] = sort(tval{i},'descend');
    plot(a,'.-'); hold on;
end
xlabel('Electrode'); ylabel('T Value');
legend(conditions);
suptitle([subject ': Max to Min T Values']);

%% Plot STRFs

% Plot STRF correlations
figure;
for j = 1:length(strfs)
    plot(mean(strfs{j}.r_values),'o-'); hold on;
end
xlabel('Electrode'); ylabel('STRF Correlation');
legend(conditions(whichouts));
suptitle([subject ': STRF Correlations by Electrode']);

% Plot each condition's STRFs on one figure
for j = 1:length(strfs)
    figure;
    subinds = numSubplots(length(bestchan));
    for i = 1:length(bestchan)
        subplot(subinds(1),subinds(2),i)
        im(squeeze(mean(strfs{j}.strf(:,:,:,i)))); axis xy
        title({['ch ' num2str(bestchan(i))], ['tval: ' num2str(tval{1}(bestchan(i)))] ['corr ' num2str(mean(strfs{j}.r_values(:,i)))]})
    end
    colormap jet
    suptitle({[subject ': STRFs'], conditions{j}});
end

% Plot STRFs for each channel with different conditions
if ploteachstrf
    for i = 1:length(bestchan)
        figure;
        for j = 1:length(strfs)
            subplot(1,length(strfs),j)
            im(squeeze(mean(strfs{j}.strf(:,:,:,i)))); axis xy
            title({conditions{j} , ['STRF Corr ' num2str(mean(strfs{j}.r_values(:,i)))]});
        end
        colormap jet;
        suptitle({['Channel ' num2str(bestchan(i))], ['tval: ' num2str(tval{1}(bestchan(i)))]});
    end
end


%% Plot average phoneme responses

% Plot average spectrograms (using the first condition)
if plotavgspec
    mean_aud = units{1}.mean_aud;
    figure;
    subinds = numSubplots(size(mean_aud,3));
    t_bef = units{1}.params.window(1);
    for phn_idx = 1:size(mean_aud,3)
        subplot(subinds(1),subinds(2),phn_idx)
        im(mean_aud(:,:,phn_idx));
        hold on; line([t_bef, t_bef]*100,[0,size(mean_aud,1)],'color','k','linestyle','--'); hold off;
        axis xy; colormap(1-gray);
        caxis manual; caxis([0 0.2*max(max(max(mean_aud)))])
        title(units{1}.units{phn_idx});
    end
    suptitle('Average Spectrogram of Each Phoneme')
end

% Plot average response to each phoneme
if plotavgeachphn
    t_bef = units{1}.params.window(1);
    for i = 1:length(units)
        mean_resp = units{i}.mean_resp;
        figure;
        for phn_idx = 1:size(mean_resp,3)
            subplot(subinds(1),subinds(2),phn_idx)
            im(mean_resp(:,:,phn_idx));
            hold on; line([t_bef, t_bef]*100,[0,size(mean_resp,1)],'color','k','linestyle','--'); hold off;
            axis xy; colormap(1-gray);
            caxis manual; caxis([0 0.2*max(max(max(mean_resp)))])
            title(units{i}.units{phn_idx});
        end
        suptitle({[subject ': Average Response to Each Phoneme'], conditions{j}});
    end
end

[mannerinds,mannerorder] = arrangePhnbyManner(units{1}.units);

% Plot average responses by channel
for j = 1:length(units)
    
    t_bef = units{1}.params.window(1);    
    
    switch whichphnrespbychan
        case {'both'}
            
            figure;
            subinds = numSubplots(length(bestchan));
            for i = 1:length(bestchan) 
                subplot(subinds(1),subinds(2),i)
                im(squeeze(units{j}.mean_resp(bestchan(i),:,mannerinds))');
                hold on; line([t_bef, t_bef]*100,[0,size(units{j}.mean_resp,1)],'color','k','linestyle','--'); hold off;
                set(gca, 'YTick', [1:1:length(units{j}.units)], 'YTickLabel', units{j}.units(mannerinds))
                xlabel('Time (10ms)');
                title({['ch ' num2str(bestchan(i))], ['tval: ' num2str(tval{1}(bestchan(i)))]})
            end
            suptitle({[subject ': Average Phoneme Reponses by Channel'], conditions{j}});
            
            figure;
            subinds = numSubplots(length(bestchan));
            for i = 1:length(bestchan)
                subplot(subinds(1),subinds(2),i)
                im(squeeze(units{j}.mean_resp_compressed(bestchan(i),:,mannerinds))');
                hold on; line([t_bef, t_bef]*100,[0,size(units{j}.mean_resp_compressed,1)],'color','k','linestyle','--'); hold off;
                set(gca, 'YTick', [1:1:length(units{j}.units)], 'YTickLabel', units{j}.units(mannerinds))
                xlabel('Time (10ms)');
                title({['ch ' num2str(bestchan(i))], ['tval: ' num2str(tval{1}(bestchan(i)))]})
            end
            suptitle({[subject ': Average Phoneme Reponses by Channel after Compression'], conditions{j}});
            
        case {'norm'}
            
            figure;
            subinds = numSubplots(length(bestchan));
            for i = 1:length(bestchan)
                subplot(subinds(1),subinds(2),i)
                im(squeeze(units{j}.mean_resp(bestchan(i),:,mannerinds))');
                hold on; line([t_bef, t_bef]*100,[0,size(units{j}.mean_resp,1)],'color','k','linestyle','--'); hold off;
                set(gca, 'YTick', [1:1:length(units{j}.units)], 'YTickLabel', units{j}.units(mannerinds))
                xlabel('Time (10ms)');
                title({['ch ' num2str(bestchan(i))], ['tval: ' num2str(tval{1}(bestchan(i)))]})
            end
            suptitle({[subject ': Average Phoneme Reponses by Channel'], conditions{j}});
            
        case {'compressed'}
            figure;
            subinds = numSubplots(length(bestchan));
            for i = 1:length(bestchan)
                subplot(subinds(1),subinds(2),i)
                im(squeeze(units{j}.mean_resp_compressed(bestchan(i),:,mannerinds))');
                hold on; line([t_bef, t_bef]*100,[0,size(units{j}.mean_resp_compressed,1)],'color','k','linestyle','--'); hold off;
                set(gca, 'YTick', [1:1:length(units{j}.units)], 'YTickLabel', units{j}.units(mannerinds))
                xlabel('Time (10ms)');
                title({['ch ' num2str(bestchan(i))], ['tval: ' num2str(tval{1}(bestchan(i)))]})
            end
            suptitle({[subject ': Average Phoneme Reponses by Channel after Compression'], conditions{j}});
    end
    
end