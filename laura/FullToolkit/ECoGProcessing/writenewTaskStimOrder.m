%% Generating new stimorder files for evnt structure
soundpath = '/Users/LauraLong/Documents/Lab/ECoG Data/Task_Sounds/'
task = 'Statistical'
tag = '*incidental*abc*';
whichblock = 2;
filenames = dir([soundpath task filesep tag]);
StimOrder = {filenames.name}'

save([soundpath task '/B' num2str(whichblock) '_stimorder.mat'],'StimOrder')