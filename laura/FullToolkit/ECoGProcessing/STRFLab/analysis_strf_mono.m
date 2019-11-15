l_r=[ones(1,32) 3*(ones(1,32)) 3*(ones(1,28)) ones(1,12) 3*(ones(1,16)) ones(1,6)];

load('./angles_64.mat');
load('./angles_59.mat');
angles=cat(1,angles_64,angles_59);
anglab = unique(angles);
angorder=[4 5 1 2 3]; %
load('out.mat')
%%
% analysis STRFresp = []; stim = [];
for i=1:123
   out(i).resp=out(i).resp(1:126,:); % only 2 subj 
end
clear resp; clear stim;
for cnt1 = 1:length(out)
    resp{cnt1} = out(cnt1).resp;
    %resp{cnt1} = cat(1,resp{cnt1},out_2(cnt1).resp);
    stim{cnt1} = out(cnt1).allauds(:,:,1);%mono
end
% remove beg and end
for i=1:length(resp)
resp{i}=resp{i}(:,125:end-50);
stim{i}=stim{i}(:,125:end-50);
end
%%

% STRF analysis for mono, predicting all angles
addpath('../')
addpath('../strflab_v1.45');

if 1 % calculate the STRFs
    ind = 1:length(angles);
    allstrf = []; fitval = []; allgof = []; presp=[]; nresp=[];
    for N=1:40%470; %size(resp{1},1)%ind2 %goodchn only     %1:size(resp{1},1) %Channels 64 for subj1 and 62 for subj 2
        disp(N);
         disp(N);
          disp(N);
         % ind_temp=[];
        ind = 1:length(angles); % prachi modification
%         for h1=2:length(anglab)
%             k=find(angles==anglab(h1));
%             ind_temp=cat(1,ind_temp,k(1:2));
%         end
        ind=[1:78 80:123]; %prachi modification
      % ind1=find(angles==0); % needed only in train0
      % ind=cat(1,ind_temp,ind1(1:3));
        tmp1=[];tmp2=[];tmp3=[];
        tolval=[0.01 0.05 0.1];
        sparseval=[ 8 16 32];
       % [strf,modelparam, presp,nresp,crr,tstind] = STRFCrossValidate2(ind,ind1,stim(ind1),resp(ind1),stim,resp,N,24,100,{0,{30,0:29},'DirectFit'},{sparseval,tolval});
        [strf,modelparam, presp,nresp,crr] = STRFCrossValidate(stim(ind),resp(ind),N,40,100,{0,{30,0:29},'DirectFit'},{sparseval,tolval});
      %  [strf,modelparam, presp,nresp,crr] = STRFCrossValidate1(ind,ind1,stim(ind),resp(ind),stim(ind1),resp(ind1),N,1,100,{0,{30,0:29},'DirectFit'},{sparseval,tolval});

       allstrf(N).strf = mean(strf,3);
        
        presp1=presp; nresp1=nresp; crr1=crr; clear crr; clear presp; clear nresp;%prachi modification
        presp(1:78)=presp1(1:78); nresp(1:78)=nresp1(1:78); crr=zeros(1,123);crr(1:78)=crr1(1:78);%
        presp(79)=cell(1,1); nresp(79)=cell(1,1);%
        presp(80:123)=presp1(79:122); %
        nresp(80:123)=nresp1(79:122);%
        crr(80:123)=crr1(79:122);
         allstrf(N).presp = presp;
        allstrf(N).nresp = nresp;
        allstrf(N).crr = crr;
        allstrf(N).ps  = [mean(modelparam.p1) mean(modelparam.p2)];
        allstrf(N).gof   = modelparam.gof;
        
        for cnt1 = 1:length(anglab)
            ind = find(angles==anglab(cnt1));
          ind=setdiff(ind,79);

%          ind_temp=[];
%         for h1=2:length(anglab)
%             k=find(angles==anglab(h1));
%             ind_temp=cat(1,ind_temp,k(1:5));
%         end
%        ind1=find(angles==0); % needed only in train0
%        ind=cat(1,ind_temp);
     
       clear nresp; clear presp;
            presp = cat(2,allstrf(N).presp{ind});
            nresp = cat(2,allstrf(N).nresp{ind});
%              presp = cat(2,(allstrf(N).presp{tstind}));
%              nresp = cat(2,allstrf(N).nresp{1:length(tstind)});

            [xData, yData] = prepareCurveData( presp, nresp );
            
            ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
            [fitresult, gof] = fit( xData, yData, ft); %, 'robust', 'bisquare' );
            fitval(1,cnt1,N) = fitresult.a;
            allgof(1,cnt1,N).gof = gof;
           
            ft = fittype( 'x+b', 'independent', 'x', 'dependent', 'y' );
            [fitresult, gof] = fit( xData, yData, ft);%, 'robust', 'bisquare' );
            fitval(2,cnt1,N) = fitresult.b;
            allgof (2,cnt1,N).gof = gof;
            
            
            ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
            [fitresult, gof] = fit( xData, yData, ft );
            fitval(3,cnt1,N) = fitresult.a;
            fitval(4,cnt1,N) = fitresult.b;
            allgof(3,cnt1,N).gof = gof;
        end
    end
else
    load STRF_mono;
end


 save(['STRF_MONO_nimascopy_silrem_lij114.mat'],'allstrf','allgof','fitval')
%%
 fitval=fitval(:,:,ind2);
l_r=l_r(ind2);
%%
[~,tmp10] = sort(l_r);

figure;

subplot(2,1,1);
imagesc(squeeze(fitval(2,[4 5 1 2 3],tmp10))); %colormap(jet)
set(gca,'yticklabel',[-90 -45 0 45 90]);
set(gca,'ytick',1:5);
colormap(jet);

ylabel('Angles'); xlabel('Electrodes'); 
subplot(2,1,2);
plot(squeeze(fitval(2,[4 5 1 2 3 ],l_r==1)),'k')
hold on
plot(squeeze(fitval(2,[4 5 1 2 3 ],l_r==3)),'r')
set(gca,'xticklabel',[-90 -45 0 45 90]);
set(gca,'xtick',1:5); ylabel('Bias');
legend('Left Brain','Right Brain')
%%
gof = [];
for cnt1 = 1:size(allgof,1)
    for cnt2 = 1:size(allgof,2)
        for cnt3 = 1:size(allgof,3)
            gof(cnt1,cnt2,cnt3) = allgof(cnt1,cnt2,cnt3).gof.adjrsquare;%sse;%rmse;%adjrsquare;
        end
    end
end
gof=gof(:,:,ind2);
a = squeeze(mean(gof,2));
figure;
subplot(2,2,1);
hold off;
bar(mean(abs(a')));
hold on;
errorbar(mean(abs(a')),std(a,[],2)/sqrt(size(a,2)),'.r');
set(gca,'xtick',1:3);
set(gca,'xticklabel',['Gain';'Bias';'Both']);
subplot(2,2,2);
scatter(a(1,:),a(3,:));
xlabel('gain modulation');
ylabel('gain and bias modulation');
line([0 0.3],[0 0.3],'linestyle','--','color','k');

subplot(2,2,3);
scatter(a(2,:),a(3,:));
xlabel('bias modulation');
ylabel('gain and bias modulation');
line([0 0.3],[0 0.3],'linestyle','--','color','k');

subplot(2,2,4);
scatter(a(1,:),a(2,:));
xlabel('gain modulation');
ylabel('bias modulation');
line([0 0.3],[0 0.3],'linestyle','--','color','k');
%% zscore
addpath('supporting files');
for i=1:size(fitval,3)
   fitval(2,:,i)=mapstd(squeeze(fitval(2,:,i)));
   
end

% divide by std
for i=1:size(fitval,3)
fitval(2,:,i)=fitval(2,:,i)./std(fitval(2,:,i));
end

% divide by abs max
for i=1:size(fitval,3)
fitval(2,:,i)=fitval(2,:,i)./abs(max(fitval(2,:,i)));
end
