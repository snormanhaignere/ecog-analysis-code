function [ecog,stim] = loadHTKDataMatrix(dpath,cond,channels,audchannel)

ecogpath = [dpath filesep cond];
stimpath = [dpath filesep 'analog'];

ecog = [];
for i = 1:length(channels)
    ecog(i,:) = readhtk([ecogpath filesep 'Ch' num2str(channels(i)) '.htk']);
end


stim = readhtk([stimpath filesep 'a' num2str(audchannel) '.htk']);



end
