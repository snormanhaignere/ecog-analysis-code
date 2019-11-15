function [strf, outmodelParam] = STRFestimate(stim,resp,fs,params,stparams)
% stim: cell array of stimulus freq * time
% resp: cell array of 1 * time
% params(1) xvalFlag:   1: perform cross validation for each stim/resp
%                       0: just estimate the strf from the whole thing.
%                       [-1 stimnumber] only cross xvalidate stimnumber
%                       2: do the jacknife
% params{2}: {1}: strfLength
%            {2}: strfDelays
%
% params{3}: 'DirectFit', 'GradDecent', 'CordDecent'
%
%   stparams(1): sparsness values
%   stparams(2): tolerance values

% wrapper for STRFlab, Nima, nimail@gmail.com

if ~exist('params','var') || isempty(params{1})
    xvalFlag = 0; % no cross validation
else
    xvalFlag = params{1};
end
if ~exist('params','var') || length(params)<2
    strfLength = 40;
    strfDelays = 0:(strfLength-1);
else
    strfLength = params{2}{1};
    strfDelays = params{2}{2};
end
if ~exist('params','var') || length(params)<3
    STRFMethod = 'DirectFit';
else
    STRFMethod = params{3};
end
if ~exist('fs','var')
    fs = 100;
end
if ~exist('stparams','var') || isempty(stparams)
    tolvals = [];
    sparsevals = [];
else
    tolvals = stparams{2};
    sparsevals = stparams{1};
end
opt = [];

%% make sure strflab functions are in matlab path, in case they aren't already
strflabDir = get_function_dir('strflab_tutorial_9_Auditory_Example');
if isempty(strflabDir)
    error('Cannot find strflab directory!');
end
addpath(genpath(strflabDir))
%%
tempPreprocDir = tempname();
[s,mess,messid] = mkdir(tempPreprocDir);
preprocStimParams.outputDir = tempPreprocDir;
%
stimInfo = struct;
stimInfo.sampleRate = fs;
stimInfo.numStimFeatures = size(stim{1},1);
%
respInfo = struct;
wholeStim = [];
wholeResponse = [];
groupIndex = [];
for cnt1 = 1:length(stim)
    wholeStim = cat(1,wholeStim,stim{cnt1}');
    wholeResponse = cat(2,wholeResponse,resp{cnt1});
    stimInfo.stimLengths(cnt1) = size(stim{cnt1},2)/fs;
    respInfo.responseLengths(cnt1) = length(resp{cnt1})/fs;
    groupIndex = cat(2,groupIndex,repmat(cnt1,[1 size(stim{cnt1},2)]));
end
%% Initialize strflab global variables with our stim and responses
global globDat
strfData(wholeStim, wholeResponse, groupIndex);
%% Initialize a linear model that extends 75ms back in time
modelParams = linInit(stimInfo.numStimFeatures, strfDelays);
%% pick training and early stopping datasets for Gradient and Coordinate Descent
N = floor(length(stim)*0.95);
trainingGroups = 1:N;
earlyStoppingGroups = N+1:length(stim);

trainingIndex = findIdx(trainingGroups, groupIndex);
earlyStoppingIndex = findIdx(earlyStoppingGroups, groupIndex);

switch STRFMethod
    case 'DirectFit'
        
        %% pick training datasets for DirectFit
        switch xvalFlag(1)
            case 0,
                trainingGroups = 1:length(stim);
                trainingIndex = findIdx(trainingGroups, groupIndex);
                %% Initialize and run the DirectFit training routine
                optOptions = trnDirectFit();
                a = optOptions.tolerances;
%                 optOptions.tolerances = a(1:4);
                optOptions.stimSampleRate = fs;
                optOptions.respSampleRate = fs;
                if ~isempty(tolvals)
                    optOptions.tolerances = tolvals;
                    optOptions.sparsenesses = sparsevals;
                end
                fprintf('\nRunning Direct Fit training...\n');
                [modelParamsDF, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);
                strf = modelParamsDF.w1;
                outmodelParam = modelParamsDF;
                outmodelParam.fs = fs;
            case 2 % do the jacknife, meaning: leave one out, perform the test, and get the mean and variance:
                trainingGroups = 1:length(stim);
                for cnt1 = trainingGroups
                    thistrain = trainingGroups;
                    thistrain(thistrain==cnt1) = [];
                    trainingIndex = findIdx(thistrain, groupIndex);
                    %% Initialize and run the DirectFit training routine
                    optOptions = trnDirectFit();
                    a = optOptions.tolerances;
                    optOptions.tolerances = a(1:4);
                    optOptions.stimSampleRate = fs;
                    optOptions.respSampleRate = fs;
                    if ~isempty(tolvals)
                        optOptions.tolerances = tolvals;
                        optOptions.sparsenesses = sparsevals;
                    end
                    fprintf('\nRunning Direct Fit training...\n');
                    [modelParamsDF, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);
                    strf(:,:,cnt1) = modelParamsDF.w1;
                    modelParamsDF.fs = fs;
                    outmodelParam(cnt1) = modelParamsDF;
                end
            case 1
                trnper = ceil(length(stim)/20); % number of examples kept for prediction
                alltrainings = 1:length(stim);
                for cnt1 = 1:length(stim)/trnper
                    trainingGroups = 1:length(stim);
                    testrange = (cnt1-1)*trnper+1:cnt1*trnper;
                    testGroups = alltrainings(testrange);
                    trainingGroups = alltrainings;
                    trainingGroups(testrange) = [];
                    trainingIndex = findIdx(trainingGroups,groupIndex);
                    testIndex = findIdx(testGroups,groupIndex);
                    % now calculate the strf:
                    optOptions = trnDirectFit();
                    optOptions.stimSampleRate = fs;
                    optOptions.respSampleRate = fs;
                    if ~isempty(tolvals)
                        optOptions.tolerances = tolvals;
                        optOptions.sparsenesses = sparsevals;
                    end
                    fprintf('\nRunning Direct Fit training...\n');
                    [modelParamsDF, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);
                    strf(:,:,cnt1)              = modelParamsDF.w1;
                    outmodelParam = modelParamsDF;
                    % now predictions:
                    [modelParamsDF, resp] = strfFwd(modelParamsDF, testIndex);
                    opt(cnt1).resp = resp;
                end
            case -1
                trainingGroups = 1:length(stim);
                trainingGroups(xvalFlag(2))=[];
                testGroups = xvalFlag(2);
                trainingIndex = findIdx(trainingGroups,groupIndex);
                testIndex = findIdx(testGroups,groupIndex);
                % now calculate the strf:
                optOptions = trnDirectFit();
                optOptions.stimSampleRate = fs;
                optOptions.respSampleRate = fs;
                if ~isempty(tolvals)
                    optOptions.tolerances = tolvals;
                    optOptions.sparsenesses = sparsevals;
                end
                fprintf('\nRunning Direct Fit training...\n');
                [modelParamsDF, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);

                strf = modelParamsDF.w1;
                outmodelParam = modelParamsDF;
                outmodelParam.fs = fs;

%                [modelParamsDF, resp] = strfFwd(modelParamsDF, testIndex);
%                 opt.resp = resp;
                
        end
    case 'GradDecent'
        
        %% Initialize and run Gradient Descent w/ early stopping
        optOptions = trnGradDesc();
        optOptions.display = 1;
        optOptions.maxIter = 300;
        optOptions.stepSize = 2e-6;
        optOptions.earlyStop = 1;
        optOptions.gradNorm = 1;
        
        fprintf('\nRunning Gradient Descent training...\n');
        [modelParams, optOptions] = strfOpt(modelParams, trainingIndex, optOptions, earlyStoppingIndex);
        strf = modelParams.w1;
        strf = squeeze(strf);
        outmodelParam = modelParams;
        %         opt.bestTol = modelParams.bestTol;
        %         opt.bestSparseness = modelParams.bestSparseness;
        %         opt.bestInfoVal = modelParams.bestInfoVal;
    case 'CordDecent'
        %% Initialize and run Coordinate Descent w/ early stopping
        optOptions = trnGradDesc();
        optOptions.display = 1;
        optOptions.maxIter = 300;
        optOptions.stepSize = 1e-4;
        optOptions.earlyStop = 1;
        optOptions.coorDesc = 1;
        
        fprintf('\nRunning Coordinate Descent training...\n');
        [modelParams, optOptions] = strfOpt(modelParams, trainingIndex, optOptions, earlyStoppingIndex);
        strf = modelParams.w1;
        outmodelParam = modelParams;
        %         opt.bestTol = modelParams.bestTol;
        %         opt.bestSparseness = modelParams.bestSparseness;
        %         opt.bestInfoVal = modelParams.bestInfoVal;
end
return;

%% split original spike trials into two PSTHs for purposes of validation
preprocRespParams.split = 1;
[wholeSplitResponse, respSplitInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);


%% compute responses to validation data for each STRF
validationGroups = [19 20];
respReal = [];
respRealHalf1 = [];
respRealHalf2 = [];
respDF = [];
respGradDesc = [];
respCoorDesc = [];
for k = 1:length(validationGroups)
    g = validationGroups(k);
    gIndx = find(globDat.groupIdx == g);
    stim = globDat.stim(gIndx);
    resp = globDat.resp(gIndx);
    respH1 = wholeSplitResponse(1, gIndx);
    respH2 = wholeSplitResponse(2, gIndx);
    
    [modelParamsDF, resp1] = strfFwd(modelParamsDF, gIndx);
    [modelParamsGradDesc, resp2] = strfFwd(modelParamsGradDesc, gIndx);
    [modelParamsCoorDesc, resp3] = strfFwd(modelParamsCoorDesc, gIndx);
    
    respReal = [respReal resp];
    respRealHalf1 = [respRealHalf1 respH1];
    respRealHalf2 = [respRealHalf2 respH2];
    respDF = [respDF resp1'];
    respGradDesc = [respGradDesc resp2'];
    respCoorDesc = [respCoorDesc resp3'];
end


%% rescale model responses
respDF = (respDF / max(respDF))*max(respReal);
respGradDesc = (respGradDesc / max(respGradDesc))*max(respReal);
respCoorDesc = (respCoorDesc / max(respCoorDesc))*max(respReal);


%% Compute performance on validation datasets for each STRF
avgNumTrials = mean(respInfo.numTrials); %taking mean isn't necessary here
infoFreqCutoff = 90; %Hz
infoWindowSize = 0.500; %500ms
[cBound, cDF] = compute_coherence_full(respDF, respReal, respRealHalf1, respRealHalf2, stimInfo.sampleRate, avgNumTrials, infoFreqCutoff, infoWindowSize);
[cBound, cGradDesc] = compute_coherence_full(respGradDesc, respReal, respRealHalf1, respRealHalf2, stimInfo.sampleRate, avgNumTrials, infoFreqCutoff, infoWindowSize);
[cBound, cCoorDesc] = compute_coherence_full(respCoorDesc, respReal, respRealHalf1, respRealHalf2, stimInfo.sampleRate, avgNumTrials, infoFreqCutoff, infoWindowSize);


%% plot STRFS
figure; hold on;

subplot(3, 1, 1); hold on;
imagesc(strfDelays, stimInfo.f, modelParamsDF.w1); axis tight;
absmax = max(max(abs(modelParamsDF.w1)));
caxis([-absmax absmax]);
colorbar;
title(sprintf('Direct Fit | bias=%f', modelParamsDF.b1));

subplot(3, 1, 2); hold on;
imagesc(strfDelays, stimInfo.f, squeeze(modelParamsGradDesc.w1)); axis tight;
absmax = max(max(abs(modelParamsGradDesc.w1)));
caxis([-absmax absmax]);
colorbar;
title(sprintf('Gradient Descent | bias=%f', modelParamsGradDesc.b1));

subplot(3, 1, 3); hold on;
imagesc(strfDelays, stimInfo.f, squeeze(modelParamsCoorDesc.w1)); axis tight;
absmax = max(max(abs(modelParamsCoorDesc.w1)));
caxis([-absmax absmax]);
colorbar;
title(sprintf('Coordinate Descent | bias=%f', modelParamsCoorDesc.b1));


%% display predictions
displayPredictions = 1;
if displayPredictions
    
    tInc = 1 / stimInfo.sampleRate;
    t = 0:tInc:(length(respReal)-1)*tInc;
    
    figure; hold on;
    
    subplot(3, 1, 1); hold on;
    plot(t, respReal, 'k-');
    plot(t, respDF, 'r-');
    axis([min(t) max(t) 0 1]);
    title('Direct Fit');
    legend('Real', 'Model');
    
    subplot(3, 1, 2); hold on;
    plot(t, respReal, 'k-');
    plot(t, respGradDesc, 'r-');
    axis([min(t) max(t) 0 1]);
    title('Gradient Descent');
    legend('Real', 'Model');
    
    subplot(3, 1, 3); hold on;
    plot(t, respReal, 'k-');
    plot(t, respCoorDesc, 'r-');
    axis([min(t) max(t) 0 1]);
    title('Coordinate Descent');
    legend('Real', 'Model');
    
end


%% plot coherences and information
figure; hold on;
plot(cBound.f, cBound.c, 'k-', 'LineWidth', 2);
plot(cDF.f, cDF.c, 'b-');
plot(cGradDesc.f, cGradDesc.c, 'g-');
plot(cCoorDesc.f, cCoorDesc.c, 'r-');
axis tight;
title('Coherences for Validation Set');
legend('Upper Bound', 'Direct Fit', 'Grad Desc', 'Coord Desc');


%% compute performance ratios
perfDF = cDF.info / cBound.info;
perfGradDesc = cGradDesc.info / cBound.info;
perfCoorDesc = cCoorDesc.info / cBound.info;

fprintf('Performance Ratios:\n');
fprintf('\tDirect Fit: %f\n', perfDF);
fprintf('\tGradient Descent: %f\n', perfGradDesc);
fprintf('\tCoordinate Descent: %f\n', perfCoorDesc);
