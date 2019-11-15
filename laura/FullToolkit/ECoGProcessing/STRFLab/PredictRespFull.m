function resp=PredictRespFull (strf,stims)
% function resp=predictresp (strf,stims)
% strf: spectro-temporal response field. Freq*Time
% stims: set of stimulus. Freq*Time


% if size(strf,2)<25, strf = resample(strf', 100, 52)';end
% for cnt1=1:size(stims,3) % number of stims
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
% alternative way:
%     m = size(stim,2) + size(strf,2);
%     a = fft(stim,m,2).*fft(strf,m,2);
%     temp = ifft(a,m,2);
% end alternative
temp(:,size(stim,2)+1:end)=[];
temp = mean(temp,1);
if exist('NumNum','var')
    for cnt1 = 1:length(NumNum)
        resp{cnt1} = temp(1:NumNum(cnt1));
        temp(1:NumNum(cnt1))=[];
    end
else
    resp = mean(temp,1);
end
temp=[];
% end