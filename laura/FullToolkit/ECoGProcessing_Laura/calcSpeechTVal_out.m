function [tval,params] = calcSpeechTVal_out(out,numframes,avgrepeats,plotvals)

disp('Running speech vs. silence significance test...');

% Params
if ~exist('numframes','var') || isempty(numframes)
    numframes = 20; dispDefaultMessage(numframes,'numframes');
end
if ~exist('plotvals','var') || isempty(plotvals)
    plotvals = 0;
end
if ~exist('avgrepeats','var') || isempty(avgrepeats)
    avgrepeats = 1; dispDefaultMessage(avgrepeats,'avgrepeats');
end

% Set/find params
bef = out(1).befaft(1); fs = out(1).dataf; % find before and dataf from out

% Find responses in both conditions
nresp = []; presp = [];
for i = 1:length(out)
    nt = 10 + floor((bef-.1)*fs*rand(1,numframes));
    pt = floor( (bef+.5)*fs + 0.5*fs*rand(1,numframes));
    if avgrepeats
        tmp = mean(out(i).resp,3);
        nresp = [nresp tmp(:,nt)];
        presp = [presp tmp(:,pt)];
    else
        for j = 1:size(out(i).resp,3)
            tmp = out(i).resp(:,:,j);
            nresp = [nresp tmp(:,nt)];
            presp = [presp tmp(:,pt)];
        end
    end
end

% Run T-test at each electrode
tval = [];
for i = 1:size(nresp,1)
    [~,~,~,tmp] = ttest2(nresp(i,:),presp(i,:));
    tval(i) = abs(tmp.tstat);
end

if plotvals
    % Display t values
    figure;
    plot(tval,'*-');
    title('Speech vs. Silence t Values by Electrode');
    ylabel('t values');
    xlabel('electrode');
    
    % Report best 10 t values
    [besttvals,bestelecs] = sort(tval,'descend');
    disp('Best ten electrodes:');
    disp(bestelecs(1:10));
end

params.numframes = numframes;
params.avgrepeats = avgrepeats;

disp('...done.');

end