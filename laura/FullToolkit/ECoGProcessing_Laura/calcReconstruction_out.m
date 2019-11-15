function [reconstruction, params] = calcReconstruction_out(out,channels,reconstructtype,rlags,avgrepeats,gap,compressresp,nC,nDS_freq,nDS_time,plotvals)
% [reconstruction] = calcReconstruction_out(out,channels,reconstructtype,rlags,stepsize,gap,compressresp,nC,nDS_freq,nDS_time,plotvals)
% Wrapper for StimuliReconstruction
%
% reconstructtype options: 
%       epochs- crossvalidate by epochs (defaults to 1/5 of the data at a time)
%       leaveoneout- crossvalidate by entry in the out structure
%       split- separate test/train sets; out is cell array {train test}
%       nocrossval- test and train on same data, no cross-validation
%
%
% Written by Laura, NAPLab April 2017

% Major params
if ~exist('channels','var') || isempty(channels)
    channels = 1:size(out(1).resp,1);
end
if ~exist('plotvals','var') || isempty(plotvals)
    plotvals = 0;
end      

% Extract params
if ~exist('avgrepeats','var') || isempty(avgrepeats)
    avgrepeats = [];
end
if ~exist('gap','var') || isempty(gap)
    gap = [];
end
if ~exist('nC','var') || isempty(nC)
    nC = [];
end
if ~exist('nDS_freq','var') || isempty(nDS_freq)
    nDS_freq = [];
end
if ~exist('nDS_time','var') || isempty(nDS_time)
    nDS_time = [];
end
if ~exist('compressresp','var') || isempty(compressresp)
    compressresp = [];
end

% Reconstruction params
if ~exist('reconstructtype','var') || isempty(reconstructtype)
    reconstructtype = 'epochs'; dispDefaultMessage(reconstructtype,'reconstructtype');
end
if ~exist('rlags','var') || isempty(rlags)
    rlags = -35:1:20; dispDefaultMessage(rlags,'rlags');
end


%% Run appropriate type of reconstruction

switch reconstructtype
        
    % Cross-validate by epoch
    case 'epochs' 
        
        disp('Running epoch-crossvalidated reconstruction...');
        
        % Extract data
        [stim, resp] = extractRespStim_out(out,channels,'array',avgrepeats,gap,compressresp,nC,nDS_freq,nDS_time);
        
        % Decide on step size
        numepochs = input('How many epochs? ');
        stepsize = ceil(size(stim,2)/numepochs);
        
        % Divide into epochs and set up vars
        epochs = [1:stepsize:size(resp,2) size(resp,2)]; % chunks
        rngs = cell(1,length(epochs)-1); % ranges
        allstim_test = zeros(size(stim)); % all test stimuli
        allrstim = zeros(size(stim)); % all reconstructed stimuli
        recons = zeros(1,length(epochs)-1); % reconstruction correlations
        g = cell(1,length(epochs)-1); % gs
        
        % Loop over epochs
        for i = 1:length(epochs)-1;
            
            % Set up train/test for this epoch
            rng = epochs(i):epochs(i+1);
            rngs{i} = rng;
            resp_test = resp(:,rng);
            resp_train = resp;
            resp_train(:,rng) = [];
            stim_test = stim(:,rng); 
            stim_train = stim;
            stim_train(:,rng) = [];
            
            % Run reconstruction
            [g{i},rstim] = StimuliReconstruction(stim_train,resp_train',resp_test',[],rlags);
            allstim_test(:,rng) = stim_test;
            allrstim(:,rng) = rstim;
            recons(i) = corrnum(stim_test,rstim);
            disp(['   corr sample ' num2str(rng(1)) ' to ' num2str(rng(end)) ': ' num2str(recons(i))]);
            
            % Set up save vars
            rstims{i} = rstim;
            stim_tests{i} = stim_test;
            
        end
        
        % Display results
        allrecon = corrnum(allstim_test,allrstim);
        disp(['Correlation for All Stimuli: ' num2str(allrecon)]);
        
        % Plot reconstruction accuracies if desired
        if plotvals
            figure;
            plot(recons,'*-');
            xlabel(['Data Chunk (' num2str(stepsize) ' Samples)']);
            ylabel('Reconstruction Coefficient');
            ylim([min(recons)*1.05 1]);
            title(['Reconstruction Accuracies by Data Chunk (' num2str(allrecon) ' for All Stimuli)']);
        end
        
      
    % Cross-validate by out structure entry (leave one out)
    case 'leaveoneout'
        
        disp('Running leave-one-out-crossvalidated reconstruction...');
        
        % Extract stim and resp as cell arrays
        [stim, resp] = extractRespStim_out(out,channels,'cell',avgrepeats,gap,compressresp,nC,nDS_freq,nDS_time); % reload as cell instead
        
        % Set up variables
        allstim_test = []; % all test stimuli
        allrstim = []; % all reconstructed stimuli
        recons = zeros(1,length(stim)); % reconstruction correlations
        g = cell(1,length(stim)); % gs
        
        % Loop over the entries
        for i = 1:length(stim)
            
            % Separate train and test
            resp_test = resp{i};
            resp_traincell = resp; resp_traincell(i) = [];
            stim_test = stim{i}; 
            stim_traincell = stim; stim_traincell(i) = [];
            resp_train = []; stim_train = [];
            for j = 1:length(resp_traincell)
                resp_train = cat(2,resp_train,resp_traincell{j});
                stim_train = cat(2,stim_train,stim_traincell{j});
            end
            
            % Run reconstruction
            [g{i},rstim] = StimuliReconstruction(stim_train,resp_train',resp_test',[],rlags);
            allstim_test = cat(2,allstim_test,stim_test);
            allrstim = cat(2,allrstim,rstim);
            recons(i) = corrnum(stim_test,rstim);
            disp(['   held out entry ' num2str(i) ': ' num2str(recons(i))]);
            
            % Set up save vars
            rstims{i} = rstim;
            stim_tests{i} = stim_test;
            
        end
        
        
        % Display results
        allrecon = corrnum(allstim_test,allrstim);
        disp(['Correlation for All Stimuli: ' num2str(allrecon)]);
        
        % Plot reconstruction accuracies if desired
        if plotvals
            figure;
            plot(recons,'*-');
            xlabel(['Out Structure Entry']);
            ylabel('Reconstruction Coefficient');
            ylim([min(recons)*1.05 1]);
            title(['Reconstruction Accuracies by Data Chunk (' num2str(allrecon) ' for All Stimuli)']);
        end
        
       
    % Different train / test sets
    case 'split'
        
        disp('Running split reconstruction...');
        
        % Get training and test sets
        [stim_train,resp_train] = extractRespStim_out(out{1},channels,'array',avgrepeats,gap,compressresp,nC,nDS_freq,nDS_time);
        [stim_test,resp_test] = extractRespStim_out(out{2},channels,'array',avgrepeats,gap,compressresp,nC,nDS_freq,nDS_time);
        
        % Run reconstruction and find correlation
        [g,rstim] = StimuliReconstruction(stim_train,resp_train',resp_test',[],rlags);
        recons = corrnum(stim_test,rstim);
        
        disp(['Correlation: ' num2str(recons)]);
    
        
    % Train/test without crossvalidation    
    case 'nocrossval' % test and train on the same data
        
        disp('Running unvalidated reconstruction...');
        
        % Extract data
        [stim, resp] = extractRespStim_out(out,channels,'array',avgrepeats,gap,compressresp,nC,nDS_freq,nDS_time);
        
        % Run reconstruction and find correlation
        [g,rstim] = StimuliReconstruction(stim,resp',resp',[],rlags);        
        recons = corrnum(stim,rstim);
        
        disp(['Correlation: ' num2str(recons)]);
        
end



%% Set up output variables

params.channels = channels;
params.reconstructtype = reconstructtype;
params.rlags = rlags;
params.avgrepeats = avgrepeats;
params.gap = gap;
params.compressresp = compressresp;
params.nC = nC;
params.nDS_freq = nDS_freq;
params.nDS_time = nDS_time;

reconstruction.g = g;
reconstruction.corr = recons;
switch lower(reconstructtype)
    case 'epochs'
        params.numepochs = numepochs;
        reconstruction.rng = rngs;
        reconstruction.stim_test = stim_tests;
        reconstruction.rstim = rstims;
        reconstruction.allrstim = allrstim;
        reconstruction.allstim_test = allstim_test;
    case 'leaveoneout'
        reconstruction.stim_test = stim_tests;
        reconstruction.rstim = rstims;
        reconstruction.allrstim = allrstim;
        reconstruction.allstim_test = allstim_test;
    case 'split'
        reconstruction.stim_test = stim_test;
        reconstruction.rstim = rstim;
    case 'nocrossval'
        reconstruction.stim_test = stim;
        reconstruction.rstim = rstim;
end

disp('...done.');

end
