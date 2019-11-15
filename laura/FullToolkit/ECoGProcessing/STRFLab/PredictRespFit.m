function [presp] = PredictRespFit (strf,stims,p)
% function resp=predictresp (strf,stims)
% strf: spectro-temporal response field. Freq*Time
% stims: set of stimulus. Freq*Time
% if resp is given, then find the corr
% param{1}: range of correlation


% if size(strf,2)<25, strf = resample(strf', 100, 52)';end
% for cnt1=1:size(stims,3) % number of stims
if ~exist('resp','var') || isempty(resp)
    resp = [];
end
if ~exist('param','var') || isempty(param)
    rng = [];
else
    rng = param{1};
end
if iscell(stims)
    % also keep the lengths:
    for cnt1 = 1:length(stims)
        NumNum(cnt1) = size(stims{cnt1},2);
    end
    stim = cat(2,stims{:});
else
    stim = stims;
end
for cnt2=1:size(stim,1)
    temp(cnt2,:)=conv(stim(cnt2,:),strf(cnt2,:));
end
% now add the fit:
temp(:,size(stim,2)+1:end)=[];
temp = sum(temp,1);
% apply the fit:
if exist('p','var') && ~isempty(p)
    temp = p(1)*temp + p(2);
%     temp =  p(1)*max(temp+p(2),0) + p(3);
end

if exist('NumNum','var')
    for cnt1 = 1:length(NumNum)
        presp{cnt1} = temp(1:NumNum(cnt1));
        temp(1:NumNum(cnt1))=[];
    end
else
    presp = mean(temp,1);
end
