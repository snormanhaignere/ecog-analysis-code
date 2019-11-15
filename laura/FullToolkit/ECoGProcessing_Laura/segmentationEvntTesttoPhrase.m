function [ phraseevnt, testoldevnt ] = segmentationEvntTesttoPhrase(subject,blocks,fixstimorevnt,keeporig )
% [ phraseevnt, testoldevnt ] = segmentationEvntTesttoPhrase(subject,blocks,fixstimorevnt,keeporig )

switch fixstimorevnt
    case 'stim'
        replaced = []; % matches out stim order and tells you which test (if any) was replaced
        testnames = {testtophrase.testwav}; % list of testnames to check for
        newStimOrder = {};
        for i = 1:length(origStimOrder)
            
            % find this stimulus name and test whether it matches one of the testnames; if yes, find tells you which one
            thisStim = origStimOrder{i};
            whichtest = find(strcmp(thisStim,testnames));
            
            % if it is a test stimulus, add the phrases for that test next; if not, just add the stimulus name
            if ~isempty(whichtest)
                thisphrasewav = testtophrase(whichtest).phrasewavs;
                newStimOrder = cat(2,newStimOrder,thisphrasewav);
                replaced = cat(2,replaced,whichtest*ones(1,length(thisphrasewav)));
            else
                newStimOrder = cat(2,newStimOrder,thisStim);
                replaced = cat(2,replaced,0);
            end
            
            
        end
        
        StimOrder = newStimOrder;
        save(['/Users/LauraLong/Documents/Lab/ECoG Data/' subject '/processed/B' num2str(block) '/Stimulus/StimOrder.mat'],'StimOrder')

        
    case 'evnt'
        
        
end



end

