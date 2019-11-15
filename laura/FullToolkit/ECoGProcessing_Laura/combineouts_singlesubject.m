%% Select which stats you're combining
subject = '062_NY668';
whichcond = 'cleaned_highgamma_norm';
savetag = '';
dpath = ['/Users/LauraLong/Documents/Lab/ECoG Data' filesep subject filesep 'processed']; % data folder
savename = 'scrambling';
savetag = '';
savepath = [dpath filesep savename];

blocklist = [1 3 4];


%% Make save directory 

if ~exist(savepath,'dir')
    mkdir(savepath);
end

%% Load all stat files

allouts = ([]); alloutfiles = cell(1,length(blocklist));
for i = 1:length(blocklist)
    
    outdir = [dpath filesep 'B' num2str(blocklist(i))];
    outfile = dir([outdir filesep '*out*' whichcond '*.mat' ]);
    if length(outfile)>1
        for j = 1:length(outfile)
            disp(outfile(j).name);
        end
        whichoutfile = input('Which out should be kept? ');
        outfile = outfile(whichoutfile);
    end
    outfile = [outdir filesep outfile.name];
    alloutfiles{i} = outfile;
    thisout = load(outfile);
    allouts{i} = thisout.out;
    
end


%% Combine structs using custom function
[out] = combineOutStructures_singlesubject(allouts,alloutfiles,blocklist);


%% Save combined data

savefilename = [savepath filesep 'out_' subject '_' savename '_' whichcond savetag '.mat'];
save(savefilename,'out','-v7.3');