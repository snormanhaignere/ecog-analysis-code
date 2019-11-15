function [normoutordata, outmean, outvar] = zscore_out( outordata, silencefile, fname )
% Z-scores an out structure or data array (channel x time)

% Determine type of data
if isequal(class(outordata),'struct')
    datatype = 'out';
elseif isequal(class(outordata),'double')
    datatype = 'data';
else
    error('datatype not recognized');
end


    
% Get resp if necessary
switch datatype 
    case 'out'
        
        if isfield(outordata,'normalization')
            disp('');
            continuequery = input('This out already has a normalization field. Continue anyway? ')
            if ~continuequery
                normoutordata = outordata; outmean = []; outvar = [];
                return
            end
        end
        
        if ~exist(silencefile,'file') || isempty(silencefile)
            %% Concatenate data
            resp = [];
            for i = 1:length(outordata)
                for j = 1:size(resp,3)
                    resp = cat(2,resp,outordata(i).resp(:,:,j));
                end
            end
        end
        
    case 'data'
        resp = outordata;
end

% Get statistics to normalize
if ~exist(silencefile,'file') || isempty(silencefile) % If no silence file, take stats from out structure
    fprintf(['...normalizing ' datatype ' by z-score from ' datatype '...']);
    %% Calculate mean and variance
    outmean = mean(resp,2);
    outvar = std(resp,[],2); % use N-1
    statfile.params = ['normalized from ' datatype ' statistics'];
    
else % Otherwise, load the silence file 
    fprintf(['...normalizing ' datatype ' by z-score from silence file...']);
    statfile = load(silencefile);
    outmean = statfile.silmean;
    outvar = statfile.silstd;
end

%% Z-score each entry in out structure or data array

switch datatype
    case 'out'
        normoutordata = outordata; % replace normout
        for i = 1:length(outordata)
            thisresp = outordata(i).resp;
            normoutordata(i).resp = (thisresp - repmat(outmean, [1, size(thisresp,2), size(thisresp,3)])) ./ repmat(outvar, [1, size(thisresp,2), size(thisresp,3)]); % subtract mean, divide by std
        end
        [normoutordata(1:end).normalization] = deal(statfile.params);
    case 'data'
        normoutordata = (resp - repmat(outmean, [1, size(resp,2), size(resp,3)])) ./ repmat(outvar, [1, size(resp,2), size(resp,3)]); % subtract mean, divide by std
end


fprintf('done.\n');

end

