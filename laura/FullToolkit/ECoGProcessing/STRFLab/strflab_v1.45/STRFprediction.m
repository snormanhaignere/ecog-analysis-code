function [presp,crr] = STRFprediction(stim,resp,inmodelParam)
% stim: cell array of stimulus freq * time
% resp: cell array of 1 * time
if ~exist('resp','var') || isempty(resp)
    resp = [];
end
fs = inmodelParam.fs;
respInfo = struct;
wholeStim = [];
wholeResponse = [];
groupIndex = [];
stimlength = [];
for cnt1 = 1:length(stim)
    stimlength(cnt1) = size(stim{cnt1},2);
    wholeStim = cat(1,wholeStim,stim{cnt1}');
    if ~isempty(resp)
        wholeResponse = cat(2,wholeResponse,resp{cnt1});
        respInfo.responseLengths(cnt1) = length(resp{cnt1})/fs;
    end
    stimInfo.stimLengths(cnt1) = size(stim{cnt1},2)/fs;
    groupIndex = cat(2,groupIndex,repmat(cnt1,[1 size(stim{cnt1},2)]));
end

testGroups = 1:length(stim);
testIndex = findIdx(testGroups,groupIndex);
testIndex = findIdx(testGroups,groupIndex);
[modelParamsDF, tmpresp] = strfFwd(inmodelParam, testIndex);
rng1 = 1; rng2 = 0; presp = [];
for cnt1 = 1:length(stim)
    rng2 = rng2 + stimlength(cnt1);
    presp{cnt1} = tmpresp(rng1:rng2);
    rng1 = rng1 + stimlength(cnt1);
    if ~isempty(resp)
        tmp1 = resp{cnt1};
        tmp2 = presp{cnt1};
        tmp2(isnan(tmp2))=0;
        crr(cnt1) = corrnum(tmp1,tmp2);
    end
end