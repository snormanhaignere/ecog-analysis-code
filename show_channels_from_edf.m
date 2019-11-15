function sig = show_channels_from_edf(exp, subjid, r, chnames)

% Plot specific channels from EDF file. Useful for example for detecting
% the trigger/audio channels.
% 
% 2019-11-11: Commented Sam NH

global root_directory

% directory for this project
project_directory = [root_directory '/' exp];

% read in specified channels
edf_file = [project_directory '/data/ECoG-EDF/' subjid '/r' num2str(r) '.edf'];
[~, sig] = edfread(edf_file, 'targetSignals', chnames);

% plot
figure;
for i = 1:length(chnames)
    subplot(length(chnames), 1, i);
    plot(sig(i,:));
    title(chnames{i});
end

