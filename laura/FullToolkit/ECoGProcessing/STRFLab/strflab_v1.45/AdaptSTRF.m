function param = AdaptSTRF(param,st,rp,rng)
if ~exist('rng','var') || isempty(rng)
    rng = 90:245;
end
h = param.w1;
h2 = h;
crr1 = [];
thismodel = param;
besth = h;
PDind = 1:length(st);
CVind = 1:length(st);
PDind = [1:10 13:22];
CVind = [11 12 23 24];
[presp,crr1] = predictresp(h2,st(CVind),rp(CVind),{rng});
initCC = (crr1);
maxCC = initCC;
disp(['initial CC: ' num2str(mean(maxCC))]);
initCC = maxCC;
subplot(2,1,1);
im(-param.w1);
subplot(2,1,2);
for cnt2 = 1:10
    %     PDind = randperm(length(st));
    %     CVind = PDind(1:ceil(length(st)*0.05));
    %     PDind(1:ceil(length(st)*0.05))=[];
   
    newcrr = []; newh = [];
    for cnt1 = 1:150
        newh(:,:,cnt1) = h2 + h2.*(randn(size(h2))-1)/2;%(1+cnt2.^.3);
%         newh(:,:,cnt1) = h2 + max2(h2)*ones(size(h2)).*(randn(size(h2))-1)/100;%(1+cnt2.^.3);
        thismodel.w1 = newh(:,:,cnt1);
        [presp,crr] = predictresp(newh(:,:,cnt1),st(PDind),rp(PDind),{rng});
        %         [presp,crr] = STRFprediction(st,rp,thismodel);
        newcrr(cnt1) = mean(crr(1:end));
    end
    newcrr(newcrr<0) =0 ;
    if ~isempty(find(isnan(newcrr))) || ~isempty(find(isnan(h2)))
        nn = 0;
    end
    if max(newcrr)>0
        newcrr = newcrr.^20;
        newcrr = newcrr - min(newcrr);
        newcrr = newcrr/std(newcrr);
        newcrr = (permute(newcrr,[1 3 2]));
        newest = newh.*repmat(newcrr,[size(newh,1) size(newh,2) 1]);
        h2 = mean(newest,3);
        %     h2 = g2/max2(g2);
        thismodel.w1 = h2;
        try
            im(-h2);drawnow;
        catch
        end
        [presp,crr1] = predictresp(h2,st(CVind),rp(CVind),{rng});
        thisCC = crr1;
        disp(['iter: ' num2str(cnt2) ' crr:' num2str(mean(thisCC)) ...
            '    ' num2str(mean(crr1(1:length(crr1)/2))) ' ' ...
            num2str(mean(crr1(length(crr1)/2+1:end)))]);
        if mean(thisCC)>mean(maxCC)
            maxCC = thisCC;
            besth = h2;
        end
    else
        disp('skipped . . .');
    end
end
param.wa = besth;
param.initCC = initCC;
param.postCC = maxCC;
