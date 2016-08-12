function f = preprocessing(subject_List)

for sub = 1:length(subject_List)
    
    D = dir([sprintf('../data/AMC%s/', subject_List{sub}), '/*.dat']);
    num_files = length(D(not([D.isdir])));
    
    for file = 1:num_files
        
        fprintf(sprintf('\n Subject %d/%d, File %d/%d...', sub, length(subject_List), file, num_files));
                
        filename.data = sprintf('../data/AMC%s/%0.2d.dat', subject_List{sub}, file);
        [signal, states, parameters] = load_bcidat(filename.data);
        signal_size = parameters.SourceCh.NumericValue;
        if strcmp(subject_List{sub}, '040') 
            signal_size = 76; % because of some inconsistencies between channel numbers in between sessions
        end
        signal = signal(:, 1:signal_size);
        sampling_rate = parameters.SamplingRate.NumericValue;
        good_channels = remove_bad_channels(signal, sampling_rate); % removed based on 60Hz noise
        nbChannels = size(signal, 2);
        bad_channels = setdiff(1:nbChannels, good_channels);
        
        bands = [8 12; 70 170];
                
        %%
                
        signal = prep_hpfilter(double(signal), 4, 0.5, parameters.SamplingRate.NumericValue);
                
        %%
        
        channel_list = good_channels;
        for k=1:ceil(nbChannels/16)
            chs = channel_list(find((channel_list>=(k-1)*16+1) & (channel_list<=k*16)));
            
            if (~isempty(chs) && length(chs)>1)
                common_average = repmat(mean(signal(:,chs),2),1,length(chs));
                signal(:,chs) = signal(:,chs)-common_average;
            end
        end
        
        %%
        
        parfor idx_band = 1:size(bands, 1)
            
            waveform_temp = prep_bpfilter(double(signal), 6, bands(idx_band, 1), bands(idx_band, 2), parameters.SamplingRate.NumericValue);
            power_temp = hilbert(single(waveform_temp));
            power_temp = abs(power_temp);
            power_temp = power_temp.^2;
            power_temp2(:, :, idx_band) = single(resample(double(power_temp), 120, parameters.SamplingRate.NumericValue));
            
        end
                
        for idx_band = 1:size(bands, 1)
            power = power_temp2(:, :, idx_band);
            save(sprintf('../features/S%s_%02d_power_%d_%d_120Hz', subject_List{sub}, file, bands(idx_band, 1), bands(idx_band, 2)), 'power')
        end
        
        stimulusCode = single(states.StimulusCode);
        stimulusCode = downsample(stimulusCode, 10);
        save(sprintf('../features/S%s_%02d_stimulusCode_120Hz', subject_List{sub}, file), 'stimulusCode')
        
        labels = parameters.Stimuli.Value;
        save(sprintf('../features/S%s_%02d_labels', subject_List{sub}, file), 'labels')
        
        clearvars -except file sub num_files subject_List
        
    end
end
end