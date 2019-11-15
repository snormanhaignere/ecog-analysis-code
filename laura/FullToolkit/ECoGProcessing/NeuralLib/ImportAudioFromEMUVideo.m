clear all;
[files,path]    = uigetfile('*.avi','select .AVI data files', 'MultiSelect','on');
mkdir(path, 'audio');

for f = 1:length(files)
    aud_file=files{f};
    aud_file(end-3:end)='.wav';
    vid = video.MultimediaFileReader([path files{f}],'AudioOutputPort',true,'VideoOutputPort',false);
    aud = video.MultimediaFileWriter([path 'audio\' aud_file],'AudioInputPort',true,'FileFormat','WAV');
    while ~isDone(vid)
    audioFrame = step(vid);
    step(aud,audioFrame);
     end
    release(aud);
    release(vid);
    f_info = dir([path files{f}]);
    data{f,1} = f_info.date;
    audio_samples = wavread([path 'audio\' aud_file]);
    data{f,2} = mean(abs(audio_samples));
end

save(files{1}(1:end-9), 'data');

figure('Position', [50 50 1240 1024]);
plot(cell2mat(data(:,2))');
title({files{1}(1:end-9) data{1,1} data{end,1}});
[Y, M, D, H, MN, S] = cellfun(@datevec,(data(:,1)));
lab = strcat(num2str(H), ':', num2str(MN));
set(gca,'XTick',1:round(size(data,1)/20):size(data,1));
set(gca,'XTickLabel',lab(1:round(size(data,1)/20):end,:));
print (gcf, '-djpeg', files{1}(1:end-9));