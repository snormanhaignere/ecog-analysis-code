function [stim, resp, params] = extractRespStim_out(out,channels,cellorarray,avgrepeats,gap,compressresp,nC,nDS_freq,nDS_time)

disp('Extracting stim and response from out...');

if ~exist('channels','var') || isempty(channels)
    channels = 1:size(out(1).resp,1); disp('defaulting to all channels');
end

if ~exist('cellorarray','var') || isempty(cellorarray)
    cellorarray = 'array'; dispDefaultMessage(cellorarray,'cellorarray');
end
if ~exist('avgrepeats','var') || isempty(avgrepeats)
    avgrepeats = 1; dispDefaultMessage(avgrepeats,'avgrepeats');
end
if ~exist('gap','var') || isempty(gap)
    gap = 0; dispDefaultMessage(gap,'gap');
end
if ~exist('nC','var') || isempty(nC)
    nC = 0.33;
    compressstim = 1; dispDefaultMessage(nC,'nC');
elseif nC == 1
    compressstim = 0;
else
    compressstim = 1;
end
if ~exist('nDS_freq','var') || isempty(nDS_freq)
    nDS_freq = 1;
    resamp_freq = 0;
elseif nC == 1
    resamp_freq = 0;
end
if ~exist('nDS_time','var') || isempty(nDS_time)
    nDS_time = 1;
    resamp_time = 0;
end
if ~exist('compressresp','var') || isempty(compressresp)
    compressresp = 1; dispDefaultMessage(compressresp,'compressresp');
end

if size(out(1).aud,3)>1
    figure;
    for i = 1:size(out(1).aud,3)
        subplot(size(out(1).aud,3),1,i)
        imagesc(out(1).aud(:,:,i).^nC);
    end
    whichaudchan = input('Which aud channel should be used? ');
else
    whichaudchan = 1;
end

switch lower(cellorarray) % Determine output format
    
    case 'array'
        
        % Set up variables
        resp = []; stim = [];
        
        % Concatenate
        for i = 1:length(out)
            if avgrepeats % if stimulus is repeated, average
                stim = cat(2,stim,out(i).aud(:,:,whichaudchan));
                resp = cat(2,resp,mean(out(i).resp,3));
                stim = cat(2,stim,zeros(size(stim,1),gap));
                resp = cat(2,resp,zeros(size(resp,1),gap));
            else % if stimulus is repeated, concatenate
                for j = 1:length(size(resp,3))
                    stim = cat(2,stim,out(i).aud(:,:,whichaudchan));
                    resp = cat(2,resp,out(i).resp(:,:,j));
                    stim = cat(2,stim,zeros(size(stim,1),gap));
                    resp = cat(2,resp,zeros(size(resp,1),gap));
                end
            end
        end
        
        % Remove unnecessary channels
        resp = resp(channels,:);
        
        % Compressions if applicable
        if compressresp
            resp = mapstd(resp);
            resp = 10*tanh(resp/10);
        end
        if compressstim
            stim = stim.^nC;
        end
        
        % Resampling if applicable
        if resamp_time
            stim = nt_dsample(stim,nDS_time);
            resp = nt_dsample(resp,nDS_time);
        end
        if resamp_freq
            stim = nt_dsample(stim,nDS_freq);
        end
        
        
    case 'cell' % do the same, but this time as a cell array (fed to calcSTRFs, keeps the out entries separated)
        
        % Set up variables
        resp = cell(1,length(out)); stim = cell(1,length(out));
        
        for i = 1:length(out)
            
            % Concatenate
            resp{i} = mean(out(i).resp,3);
            stim{i} = out(i).aud(:,:,whichaudchan);
            
            % Compressions if applicable
            if compressresp
                resp{i} = mapstd(resp{i});
                resp{i} = 10*tanh(resp{i}/10);
            end
            if compressstim
                stim{i} = stim{i}.^nC;
            end
            
            % Remove unnecessary channels
            resp{i} = resp{i}(channels,:);
            
            % Resampling if applicable
            if resamp_time
                stim{i} = nt_dsample(stim{i},nDS_time);
                resp{i} = nt_dsample(resp{i},nDS_time);
            end
            if resamp_freq
                stim{i} = nt_dsample(stim{i},nDS_freq);
            end
            
        end
        
end

params.channels = channels;
params.avgrepeats = avgrepeats;
params.gap = gap;
params.compressresp = compressresp;
params.nC = nC;
params.nDS_freq = nDS_freq;
params.nDS_time = nDS_time;

disp('...done.');

end