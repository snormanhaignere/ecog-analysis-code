function [strf, outmodelParam, presp, nresp, crr] = STRFCrossValidate(stim,resp,chn,N,fs,params,stparams)
% in this case, stim and resp should be cell arrays!

ind = randperm(length(stim));
presp = cell(0);
Nn = ceil(length(stim)/N);
nresp = [];
p1 = [];
p2 = []; p3=[];
for cnt1 = 1:length(resp)
    nresp{cnt1} = resp{cnt1}(chn,:);
end
for cnt1 = 1:Nn
    disp([num2str(cnt1) ' out of ' num2str(Nn)]);
    trnind = ind;
    tstind = trnind((cnt1-1)*N+1:min(length(stim),cnt1*N));
    trnind((cnt1-1)*N+1:min(length(stim),cnt1*N))=[];
    stimtrn = cat(2,stim{trnind});
    resptrn = cat(2,nresp{trnind});
    stimtst = cat(2,stim{tstind});
    resptst = cat(2,nresp{tstind});
    
    [strf(:,:,cnt1),outmodelParam] = STRFestimate({stimtrn},{resptrn},fs,params,stparams);
    % now fit a&b, where r = a SRTF  + b
    tmp = PredictRespFit(strf(:,:,cnt1),stim(trnind));
    tmp1 = cat(2,tmp{:});
    tmp2 = cat(2,resp{trnind});
    tmp2 = tmp2(chn,:);
    % now find the fit:
    [xData, yData] = prepareCurveData( tmp1, tmp2 );
    
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft);
    p1(end+1) = fitresult.p1;
    p2(end+1) = fitresult.p2;
%     p3(end+1) = fitresult.c;
    allgof(cnt1) = gof;
    % and predict the test:
    tmp = PredictRespFit(strf(:,:,cnt1),stim(tstind),[fitresult.p1 fitresult.p2]);% fitresult.c]);
    presp(tstind) = tmp;
end
outmodelParam.p1 = p1;
outmodelParam.p2 = p2;
% outmodelParam.p3 = p3;
outmodelParam.gof = allgof;
crr = [];
for cnt1 = 1:length(presp)
    crr(cnt1) = corrnum(presp{cnt1},nresp{cnt1});
end