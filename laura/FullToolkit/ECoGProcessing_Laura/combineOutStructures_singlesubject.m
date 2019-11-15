function [ allout ] = combineOutStructures_singlesubject(outs,outlocs,outnames)
% [ allout ] = combineOutStructures_Laura(outs,outlocs)
%   can load either with outs being the structures or with outlocs being the
%    location of the out structures
%   outnames is the name that will be loaded into allout that helps
%   distinguish which elecs come from which

% Determine whether outs need to be loaded
if exist('outs','var') && ~isempty(outs)
    loadouts = 0;
elseif exist('outlocs','var') && ~isempty(outlocs)
    loadouts = 1;
else
    error('need outs or outlocs');
end

% Initialize with first out and grab fields
if loadouts
    allout = load(outlocs{1});
else
    allout = outs{1};
end
allfields = fields(allout);
alloutnames = {allout.name}';
numtrials = length(allout);

% % Add subject,whichout,whichoutfile fields to track which ch are from where
% numchan(1) = size(allout(1).resp,1);
% for i = 1:length(allout)
%     [allout(i).subject{1:numchan}] = deal(outnames{1});
%     allout(i).numchaneachout(1) = numchan;
%     [allout(i).whichout(1:numchan)] = deal(1);
%     if exist('outlocs','var') && ~isempty(outlocs)
%         [allout(i).whichoutfile{1:numchan}] = deal(outlocs{1});
%     end
% end


% Now add the others
for i = 2:length(outs)
    
    % Load out
    if loadouts
        thisout = load(outlocs{i});
    else
        thisout = outs{i};
    end
    % Compare fields with allfields; if it doesn't match, error
    if ~isequal(allfields,fields(thisout))
        disp(allfields); disp(fields(thisout));
        error('fields do not match, see above');
    end
    
    % Now loop over each out entry and add the new data
    for j = 1:length(thisout)
        
        disp(numtrials);
        
        
        thisname = thisout(j).name;
        nameloc = find(strcmpi(thisname,alloutnames));
        if isempty(nameloc)
            nameloc = length(allout)+1;
            allout(nameloc) = thisout(j);
            allout(nameloc).trial = allout(nameloc).trial + numtrials;
        else
            % Loop over original fields
            for k = 1:length(allfields)
                thisfield = allfields{k};
                switch thisfield
                    case 'resp'
                        allout(nameloc).(thisfield) = cat(3,allout(nameloc).(thisfield),thisout(j).(thisfield));
                    case 'audrecord'
                        allout(nameloc).(thisfield) = cat(2,allout(nameloc).(thisfield),thisout(j).(thisfield));
                    case 'chnames'
                        allout(nameloc).(thisfield) = cat(2,allout(nameloc).(thisfield),thisout(j).(thisfield));
                    case 'params'
                        allout(nameloc).(thisfield) = cat(1,allout(nameloc).(thisfield),thisout(j).(thisfield));
                    case 'trial'
                        allout(nameloc).(thisfield) = cat(2,allout(nameloc).(thisfield),thisout(j).(thisfield)+numtrials); % for trials, add the number of previous fields to make sure it's corrected for multiple blocks
                    otherwise % for others, check that they're identical therefore combinable
                        if ~isequal(allout(nameloc).(thisfield),thisout(j).(thisfield));
                            try
                                discrepancy = abs(allout(nameloc).(thisfield) - thisout(j).(thisfield));
                                if discrepancy > 1e-5
                                    error([thisfield ' discrepancy is ' num2str(discrepancy)]);
                                end
                            catch
                                disp(allout(nameloc).(thisfield)); disp(thisout(j).(thisfield));
                                error(['content of ' thisfield ' field does not match']);
                            end
                        end
                end
            end
            
        end
    end
    
    
    alloutnames = {allout.name}';
    numtrials = numtrials + length(thisout);

end

allout = sortstruct(allout,'name');

end

