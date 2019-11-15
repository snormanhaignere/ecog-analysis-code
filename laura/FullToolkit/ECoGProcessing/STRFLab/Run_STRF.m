function out = Run_STRF(stim,resp,dir,lags,nDS_freq, nDS_time,nC,k,method)
% Function for calculating STRFs or Reconstruction filters
%
% Data are automatically split into training and testing sets.
% Filters are fit using Ridge Regression. 
% Optimal ridge parameters are automatically selected to prevent overfitting. 
%
% INPUTS
% ------------------------------------------------------------------------------------------------------------
% Required
% stim: Stimulus spectrogram 
% resp: Neural response data 
%
% Optional
% dir: Direction of regression. 1 = forward (for STRF; default). -1 = backward (for Stimulus-Reconstruction).
%
% lags: Time-lags in samples
% (e.g., if fs = 100Hz, and you want time-lags spanning -100 to 400ms, lags = -10:40;
% (N.B.!!!!!!!! If you are generating reconstruction filters, the time-lags
% need to be reversed. So -100 to 400 becomes, lags = -40:10;
%
% nDS_freq: Downsample factor for spectrograms in frequency dimension (e.g., nDS = 8: 128 to 16 frequency bands)
% nDS_time: Downsample factor in time
% nC: Compression factor for spectrograms (Default = 0.33)
% k: value for k-fold cross-validation
% method: Type of regression performed. 'ridge' (Default) or 'svd'. 

% OUTPUTS 
% ------------------------------------------------------------------------------------------------------------
% strf: 
% if dir = 1: spectrotemporal receptive field (STRF)
% if dir = -1: Reconstruction filter
% (unwrapped for plotting purposes, and chosen via cross-validation)
%
% r_values: 
% if dir = 1: correlation (Pearson's r) of predicted and actual responses for each electrode 
% if dir = -1: correlation (Pearson's r) of reconstructed and actual spectrograms for each frequency band
% (all validation sets, best tolerance)
%
% strf2: 
% corresponding STRF (decoder) used to convolve with
% spectrogram (neural data) to predict neural responses (reconstruct spectrogram)
%
% prediction: 
% if dir = 1: predicted data for each cross-validation
% if dir = -1: reconstructed spectrograms for each cross-validation
%
% actual: 
% if dir = 1: actual neural response for each cross-validation
% if dir = -1: actual spectrogram for each cross-validation
%
% tol: selected tolerance value
%
% Written by James O'Sullivan
% Lab of Nima Mesgarani
% 2016

if ~exist('dir','var') | isempty(dir)
    dir = 1;
end

if ~exist('lags','var') | isempty(lags)
    lags = -10:40;
end

if ~exist('nDS_freq','var') | isempty(nDS_freq)
    nDS_freq = 8;
end
if ~exist('nDS_time','var') | isempty(nDS_time)
    nDS_time = 1;
end
if ~exist('nC','var') | isempty(nC)
    nC = 0.33;
end

if ~exist('k','var') | isempty(k)
    k = 5;
end

if ~exist('method','var') | isempty(k)
    method = 'ridge';
end

% Set up ridge parameters (lambda)
lambda = logspace(0,log10(1e5),21);

% Check dimensions
if size(stim,2) > size(stim,1)
    stim = stim';
end
if size(resp,2) > size(resp,1)
    resp = resp';
end

if size(resp,1) ~= size(stim,1)
    error('Size Mismatch')
end

% Make non-negative
stim = stim-min2(stim);

% Compress Spectrograms
stim = stim.^nC;

% Normalise
stim = (stim - mean2(stim))./std2(stim);

% LL 9/17- commented this out from James's input
% if dir == 1 
%     % Forward model. Each electrode is independent, so normalize each individually. 
%     resp = zscore(resp);
% else
%     % Backward model. Electrode are not independent, so normalize all simultaneously.
%     resp = (resp - mean2(resp))./std2(resp);
% end

% Downsample Spectrogram
stim = nt_dsample(stim',nDS_freq)';

% Downsample Time
stim = nt_dsample(stim,nDS_time);
resp = nt_dsample(resp,nDS_time);

% Number of parameters that will be fit
nParams = length(lags)*size(stim,2)+1;

% Get training and test length
clear r rM w W
TrainL = floor(size(stim,1)/k); % Train Length 
TestL = size(stim,1) - TrainL; % Test Length 

% Compare training length with number of parameters
if nParams < TrainL*10
%     warning('Too many parameters to fit given data size. Suggest reducing number of lags/frequency bands/k')
end

nFreqs = size(stim,2);
nElecs = size(resp,2);

% initialize matrices
if dir == 1
    strf = zeros(k,nFreqs,length(lags),size(resp,2));    
    strf_all = zeros(k,length(lambda),nFreqs,length(lags),size(resp,2));
    strf2 = zeros(k,nFreqs*length(lags)+1,size(resp,2));
    strf2_all = zeros(k,length(lambda),nFreqs*length(lags)+1,size(resp,2));
    r_values_all = zeros(k,length(lambda),nElecs);
    mse_all = zeros(k,length(lambda),nElecs);
    testPred = zeros(k,length(lambda),size(resp,2),TestL);
    prediction = zeros(k,size(resp,2),TestL);
else
    strf = zeros(k,nElecs,length(lags),size(stim,2));
    strf_all = zeros(k,length(lambda),nElecs,length(lags),size(stim,2));
    strf2 = zeros(k,nElecs*length(lags)+1,size(stim,2));
    strf2_all = zeros(k,length(lambda),nElecs*length(lags)+1,size(stim,2));
    r_values_all = zeros(k,length(lambda),nFreqs);
    mse_all = zeros(k,length(lambda),nFreqs);
    testPred = zeros(k,length(lambda),size(stim,2),TestL);
    prediction = zeros(k,size(stim,2),TestL);
end

testStim = zeros(k,TestL,size(stim,2));
testResp = zeros(k,TestL,size(resp,2));

% k-fold cross-validation
for i = 1:k    
    % Get training samples
    TrainSamples = (i-1)*TrainL+1:i*TrainL;
    
    % Get test samples    
    TestSamples = setdiff(1:size(stim,1),TrainSamples);
    
    % Get data from test set
    testStim(i,:,:) = stim(TestSamples,:);
    testResp(i,:,:) = resp(TestSamples,:);     
    
    % Calculate STRFs on training data
    if dir == 1
        strf2_all(i,:,:,:) = STRF(stim(TrainSamples,:),resp(TrainSamples,:),lags,lambda,method);
        lag = lagmatrix(squeeze(testStim(i,:,:)),lags);
    else
        strf2_all(i,:,:,:) = STRF(resp(TrainSamples,:),stim(TrainSamples,:),lags,lambda,method);
        lag = lagmatrix(squeeze(testResp(i,:,:)),lags);
    end        

    lag = [ones(size(lag,1),1) lag]; %#ok<AGROW>

    % Loop through each lambda value and predict data
    for j = 1:length(lambda)
        
        % Predict data via convolution with STRF
        testPred(i,j,:,:) = squeeze(strf2_all(i,j,:,:))'*lag';
        
        % Get correlation
        if dir == 1
            r_values_all(i,j,:) = diag(corr(squeeze(testPred(i,j,:,:))',squeeze(testResp(i,:,:))));
            
            % Get mean squared error
            for h = 1:nElecs
                mse_all(i,j,h) = immse(squeeze(testPred(i,j,h,:)),squeeze(testResp(i,:,h))');
            end            
            % Unwrap STRF for plotting purposes
            for f = 1:nFreqs
                strf_all(i,j,f,:,:) = strf2_all(i,j,f+1:nFreqs:end,:);
            end
        else
            r_values_all(i,j,:) = diag(corr(squeeze(testPred(i,j,:,:))',squeeze(testStim(i,:,:))));
            for h = 1:nFreqs
                mse_all(i,j,h) = immse(squeeze(testPred(i,j,h,:)),squeeze(testStim(i,:,h))');
            end            
            for f = 1:nElecs
                strf_all(i,j,f,:,:) = strf2_all(i,j,f+1:nElecs:end,:);
            end
        end            
    end
   
end

% Get STRF for best tolerance
clear r_values tol strf strf2 prediction
if dir == 1
    x = nElecs;
else
    x = nFreqs;
end

r_values = zeros(k,x);
tol = zeros(1,x);
mse = zeros(1,x);

for i = 1:x
    [mse(i),tol(i)] = min(mean(mse_all(:,:,i),1));
    r_values(:,i) = squeeze(r_values_all(:,tol(i),i));
    strf(:,:,:,i) = squeeze(strf_all(:,tol(i),:,:,i));
    strf2(:,:,i) = squeeze(strf2_all(:,tol(i),:,i));
    prediction(:,i,:) = squeeze(testPred(:,tol(i),i,:));
end

if dir == 1
    actual = permute(testResp,[1 3 2]);
else
    actual = permute(testStim,[1 3 2]);
end

% Assign variables to output structure
out.strf = strf;
out.r_values = r_values;
out.mse = mse;
out.strf2 = strf2;
out.prediction = prediction;
out.actual = actual;
out.tol_mse = tol;

% Other variables that can be saved if wanted
% out.strf_all = strf_all;
% out.strf2_all = strf2_all;
% out.r_values_all = r_values_all;
% out.mse_all = mse_all;
