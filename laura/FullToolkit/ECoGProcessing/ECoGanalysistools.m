%% Select Portions of Code to Run
avgresponse = 0;
extract = 1;
sigtest = 1;
reconstruct_crossval = 0;
reconstruct_allvals = 0;
reconstruct = 0;
strf = 1;


%% Plot average response to stimulus from the 

if avgresponse
   
    range = 1:160;
    
    data = zeros(size(out(1).resp,1),length(range),length(out));
    for i = 1:length(out)
        data(:,:,i) = mean(out(i).resp(:,range,:),3);
    end
    
    figure; 
    imagesc(mean(data,3));
    title('Average Response to Stimulus');
    ylabel('Electrodes'); xlabel('Time (10ms)');
    
    
end


%% Extract the Stim and Resp

if extract
    disp('Extracting Stim and Resp');
    
    % Set up variables
    resp = []; stim = [];
    gp = 0; % set gap between trials if desired
    
    % Loop to load stim and resp into
    for cnt1 = 1:length(out)
        disp(cnt1);
        stim = cat(2,stim,out(cnt1).aud);
        resp = cat(2,resp,mean(out(cnt1).resp,3));
        stim = cat(2,stim,zeros(size(stim,1),gp));
        resp = cat(2,resp,zeros(size(resp,1),gp));
    end
    
    % Compress resp and stim
    resp = mapstd(resp);
    resp = 10*tanh(resp/10);
    stim = stim.^.25;
    
end

%% Significance Test
% Based on assumption that response to clean and sound should have different magnitudue

if sigtest
    disp('Significance Test');
    
    % Set/find params
    NF = 20; % number of frames from response and clean
    bef = out(1).befaft(1); fs = out(1).dataf; % find before and dataf from out
    
    % Find responses in both conditions
    nresp = []; presp = [];
    for cnt1 = 1:length(out)
        disp(cnt1);
        nt = 10 + floor((bef-.1)*fs*rand(1,NF));
        pt = floor( (bef+.5)*fs + 0.5*fs*rand(1,NF));
        tmp = mean(out(cnt1).resp,3);
        nresp = [nresp tmp(:,nt)];
        presp = [presp tmp(:,pt)];
    end
    
    % Run T-test at each electrode
    tval = [];
    for cnt1 = 1:size(nresp,1)
        [~,~,~,tmp] = ttest2(nresp(cnt1,:),presp(cnt1,:));
        tval(cnt1) = abs(tmp.tstat);
    end
    
    % Display t values
    figure;
    plot(tval,'*-');
    title('t Values by Electrode');
    ylabel('t values');
    xlabel('electrode');
    
    % Report best 10 t values
    [besttvals,bestelecs] = sort(tval,'descend');
    disp('Best ten electrodes:');
    disp(bestelecs(1:10));
    
end


%% Reconstruction
if reconstruct_crossval
    
    rng = 10000:16000; % set sample range to plot (overwrites rng for waveform)
    
    resp_test = resp(:,rng);
    resp_train = resp;
    resp_train(:,rng) = [];
    
    stim_test = stim(:,rng);
    stim_train = stim;
    stim_train(:,rng) = [];
    
    % Run Reconstruction
    rlags = -35:1:20;
    [g,rstim] = StimuliReconstruction(stim_train,resp_train',resp_test',[],rlags);
    recon = corrnum(stim_test,rstim)
    
    % Plot example spectrograms for comparison
    fs = out.dataf;
    timelabels = rng/fs;
    freqlabels = 1:128;
    
    figure;
    subplot(2,1,1);
    imagesc(timelabels,freqlabels,stim_test); axis xy
    title('Original');
    subplot(2,1,2);
    imagesc(timelabels,freqlabels,rstim);
    xlabel('Time (s)'); ylabel('Frequency Band');
    title(['Reconstruction     Correlation: ' num2str(recon)]); axis xy
    xlabel('Time (s)'); ylabel('Frequency Band');
    colormap(jet);
end


if reconstruct_allvals
    
    stepsize = 5000; % epoch step size
    
    epochs = [1:stepsize:size(resp,2) size(resp,2)]; % chunks
    rngs = cell(1,length(epochs-1)); % ranges
    allstim_test = zeros(size(stim)); % all test stimuli
    allrstim = zeros(size(stim)); % all reconstructed stimuli
    recons = zeros(1,length(epochs)-1); % reconstruction correlations
    g = cell(1,length(epochs)-1); % gs
    for i = 1:length(epochs)-1;
        rng = epochs(i):epochs(i+1);
        rngs{i} = rng;
        
        resp_test = resp(:,rng);
        resp_train = resp;
        resp_train(:,rng) = [];
        
        stim_test = stim(:,rng);
        stim_train = stim;
        stim_train(:,rng) = [];
        
        % Run Reconstruction
        rlags = -35:1:20;
        [g{i},rstim] = StimuliReconstruction(stim_train,resp_train',resp_test',[],rlags);
        allstim_test(:,rng) = stim_test;
        allrstim(:,rng) = rstim;
        recons(i) = corrnum(stim_test,rstim);
        disp(['   corr sample ' num2str(rng(1)) ' to ' num2str(rng(end)) ': ' num2str(recons(i))]);
        
    end
    
    allrecon = corrnum(allstim_test,allrstim);
    disp(['Reconstruction Correlation for All Stimuli: ' num2str(allrecon)]);
    
    % Plot reconstruction accuracies
    figure;
    plot(recons);
    xlabel(['Data Chunk (' num2str(stepsize) ' Samples)']);
    ylabel('Reconstruction Coefficient');
    title(['Reconstruction Accuracies by Data Chunk (' num2str(allrecon) ' for All Stimuli)']);
    
    % Plot one segment
    whichrng = [1 12]; % pick which chunk to plot
    fs = out.dataf;
    timelabels = rng/fs;
    freqlabels = 1:128;
    
    for i = 1:length(whichrng);
        figure;
        subplot(2,1,1);
        imagesc(timelabels,freqlabels,allstim_test(:,rngs{whichrng(i)})); axis xy
        title('Original');
        subplot(2,1,2);
        imagesc(timelabels,freqlabels,allrstim(:,rngs{whichrng(i)}));
        xlabel('Time (s)'); ylabel('Frequency Band');
        title(['Reconstruction     Correlation: ' num2str(recons(whichrng(i)))]); axis xy
        xlabel('Time (s)'); ylabel('Frequency Band');
        colormap(jet);
    end
    
end


%% Reconstruction

if reconstruct
    disp('Stimulus Reconstruction');
    
    
    
    % Run Reconstruction!!
    rlags = -35:1:20;
    [g,rstim] = StimuliReconstruction(stim,resp',resp',[],rlags);
    
    
    % Find correlation!
    corrnum(stim,rstim)
    
    % Correlation per frequency
    crng = find(mean(stim,1)>0);
    corrperfreq(stim(:,crng),rstim(:,crng))
    
    
    % Convert portion of rstim to waveform
    rng = 1:500; % set sample range to turn into waveform
    
    tmp = max(rstim(:,rng),0); % half-wave rectification to ensure spectrograms are positive
    tmp = tmp.^(1/.25); % inverting the compression
    rw = aud2wav(tmp',[],[10 10 -2 0 40 1 0]); % rw is waveform of reconstructed stimulus
    
    
    % Plot example spectrograms for comparison
    rng = 1:1000; % set sample range to plot (overwrites rng for waveform)
    
    figure;
    subplot(2,1,1);
    imagesc(rstim(:,rng));
    title('Reconstruction'); axis xy
    subplot(2,1,2)
    imagesc(stim(:,rng)); axis xy
    title('Original');
    suptitle([num2str(rng(1)/100) ' to ' num2str(rng(end)/100) ' s']);
    
end


%% STRF ecog:

if strf
    disp('Finding STRFs');
    
    
    strfcorr = zeros(1,length(tval));
    
    resampstim = resample(stim,1,5);
    allstrfs = zeros(size(resampstim,1),30,length(tval));
    
    figure;
    
    for i = 1:length(tval);
        N = bestelecs(i); % select electrode
        testrng = 1:1000; % set sample range to plot for comparison
        trainrng = 1:size(stim,2);
        trainrng(testrng) = [];
        
        % Run STRF code
        disp(i);
        tolval=[.005 .01] ; sparseval=[4 6 8];
        [strfN,tmp2] = STRFestimate({resampstim(:,trainrng)},{resp(N,trainrng)},100,{0,{30,0:29},'DirectFit'},{sparseval,tolval});
        allstrfs(:,:,i) = strfN;
        allstrfparams(i) = tmp2;
        
        % Predict the response to a sound given the STRF
        presp  = predictresp(strfN,resampstim); % find predicted response given STRF and stimulus
        strfcorr(i) = corr(presp(testrng)',resp(N,testrng)');
        
        % Draw STRF
        timelabels = [1:300];
        freqlabels = [1:size(strfN,1)];
        mm = max(abs(strfN(:)));
        
        scrollsubplot(3,3,i);
        colormap jet;
        imagesc(timelabels,freqlabels,strfN,[-mm mm]);
        xlabel('Time (ms)'); ylabel('Frequency Band');
        if exist('chnames','var')
            title({[num2str(N) ': ' chnames{N} ] ; ['tval: ' num2str(tval(N)) ' corr: ' num2str(strfcorr(i))] });
        else
            title({['elec ' num2str(N)] ; ['tval: ' num2str(tval(N)) ' corr: ' num2str(strfcorr(i))] });
        end
        drawnow;
        
    end
    
    suptitle('STRFs by Electrode');
    
    strfdone = 1;
    
end

%%

% if strf
%     disp('Finding STRF');
%
%     N = 12; % select electrode
%     rng = 1:1000; % set sample range to plot for comparison
%
%     % Run STRF code
%     tolval=[.0005 .001 .005 .01] ; sparseval=[4 6 8];
%     [strfN,tmp2] = STRFestimate({stim(:,1000:end)},{resp(N,1000:end)},100,{0,{30,0:29},'DirectFit'},{sparseval,tolval});
%
%
%     % Predict the response to a sound given the STRF
%     presp  = predictresp(strfN,stim); % find predicted response given STRF and stimulus
%     strfcorr = corr(presp',resp(N,:)')
%     figure; % plot responses
%     plotm1(resp(N,rng),'b');
%     hold on;
%     plotm1(presp(rng),'r');
%     title({['Electrode ' num2str(N)] ; ['t Value: ' num2str(tval(N)) '  Correlation: ' num2str(strfcorr)] ; 'Actual response: blue           Predicted response: red'});
%
%     % Draw STRF
%     figure;
%     colormap jet;
%     im(strfN);drawnow;
%     title(['Best STRF for Electrode ' num2str(N)]);
%     title({num2str(N) ; ['tval: ' num2str(tval(N)) '           corr: ' num2str(strfcorr)] });
%
%
% end

