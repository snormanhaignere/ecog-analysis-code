%% Load data


load(['out_025_LIJ109_B85_CMRamp_highgamma.mat']); % the one with 64 trials
out2=load(['out_025_LIJ109_B86_CMRamp_highgamma.mat']);
out2=out2.out;

out=cat(1,out(:),out2(:));

load('angles_64.mat');
load('angles_59.mat');


angles=cat(1,angles_64,angles_59);

%%


for cnt1 = 1:length(out)
    %     tmp = out(cnt1).sound(:,1);
    %     tmp = resample(tmp,16000,out(1).soundf);
    %     tmp = wav2aud(tmp,[10 10 -2 0])';
    %     tmp = resample(tmp.^.25,1,8);
    %
    %     out(cnt1).newaudleft=tmp;
    
    tmp=[];
    tmp = out(cnt1).sound(:,1);
    tmp = resample(tmp,16000,out(1).soundf);
    tmp = wav2aud(tmp,[10 10 -2 0])';
    tmp = resample(tmp.^.25,1,8);
    while size((out(cnt1).resp),2) < size(tmp,2)
        tmp(:,size(tmp,2))=[];
    end
    
    out(cnt1).newaudleft=tmp;
    
    % loadload; close;
end


%%

for i=1:123
   out(i).newaudcat=cat(1,out(i).newaudleft,out(i).newaudright); 
    
end
%%

resp = [];
stim = [];

resp_c1 = [];
stim_c1 = [];
resp_c2 = [];
stim_c2 = [];
resp_c3 = [];
stim_c3 = [];
resp_c4 = [];
stim_c4 = [];
resp_c5 = [];
stim_c5 = [];


gp = 0;

% %Labels
% for i=1:64
% angles(i)=out(i).Angles;
% end


%
for cnt1 = 1:length(out)
    if angles(cnt1)==0
        stim_c1 = cat(2,stim_c1,squeeze(out(cnt1).newaudcat(:,11:end-10)));
        resp_c1 = cat(2,resp_c1,mean(out(cnt1).resp(:,11:end-10,:),3));
        stim_c1 = cat(2,stim_c1,zeros(size(stim_c1,1),gp));
        resp_c1 = cat(2,resp_c1,zeros(size(resp_c1,1),gp));
    end
    if angles(cnt1)==45
        stim_c2 = cat(2,stim_c2,squeeze(out(cnt1).newaudcat(:,11:end-10)));
        resp_c2 = cat(2,resp_c2,mean(out(cnt1).resp(:,11:end-10,:),3));
        stim_c2 = cat(2,stim_c2,zeros(size(stim_c2,1),gp));
        resp_c2 = cat(2,resp_c2,zeros(size(resp_c2,1),gp));
    end
    if angles(cnt1)==90
        stim_c3 = cat(2,stim_c3,squeeze(out(cnt1).newaudcat(:,11:end-10)));
        resp_c3 = cat(2,resp_c3,mean(out(cnt1).resp(:,11:end-10,:),3));
        stim_c3 = cat(2,stim_c3,zeros(size(stim_c3,1),gp));
        resp_c3 = cat(2,resp_c3,zeros(size(resp_c3,1),gp));
    end
    if angles(cnt1)==270
        stim_c4 = cat(2,stim_c4,squeeze(out(cnt1).newaudcat(:,11:end-10)));
        resp_c4 = cat(2,resp_c4,mean(out(cnt1).resp(:,11:end-10,:),3));
        stim_c4 = cat(2,stim_c4,zeros(size(stim_c4,1),gp));
        resp_c4 = cat(2,resp_c4,zeros(size(resp_c4,1),gp));
    end
    if angles(cnt1)==315
        stim_c5 = cat(2,stim_c5,squeeze(out(cnt1).newaudcat(:,11:end-10)));
        resp_c5 = cat(2,resp_c5,mean(out(cnt1).resp(:,11:end-10,:),3));
        stim_c5 = cat(2,stim_c5,zeros(size(stim_c5,1),gp));
        resp_c5 = cat(2,resp_c5,zeros(size(resp_c5,1),gp));
    end
    
    stim = cat(2,stim,out(cnt1).newaudcat(:,11:end-10));
    resp = cat(2,resp,mean(out(cnt1).resp(:,11:end-10,:),3));
    stim = cat(2,stim,zeros(size(stim,1),gp));
    resp = cat(2,resp,zeros(size(resp,1),gp));
    
end



%Remove blank channels

resp=resp(1:ch,:);
resp_c1=resp_c1(1:ch,:);
resp_c2=resp_c2(1:ch,:);
resp_c3=resp_c3(1:ch,:);
resp_c4=resp_c4(1:ch,:);
resp_c5=resp_c5(1:ch,:);


%%


%% STRF ecog:

for N=1:64
    tmp1=[];
    tmp2=[];
    %tolval=[.0005 .001 .005 .01] ;
    tolval=0.1;
    sparseval=[ 4 ];
    %for ang=1:5
    ang=1;
    %change angle   c1..5
    % leave out 2 trials 11:1358
    [tmp1(ang,:,:),tmp2] = STRFestimate({stim_c1(:,1359:end)},{resp_c1(N,1359:end)},100,{0,{30,0:29},'DirectFit'},{sparseval,tolval});
    
    
    
    ang=2;
    % leave out 1:780
    [tmp1(ang,:,:),tmp2] = STRFestimate({stim_c2(:,781:end)},{resp_c2(N,781:end)},100,{0,{30,0:29},'DirectFit'},{sparseval,tolval});
    
    ang=3;
    % leave out 1:963
    [tmp1(ang,:,:),tmp2] = STRFestimate({stim_c3(:,964:end)},{resp_c3(N,964:end)},100,{0,{30,0:29},'DirectFit'},{sparseval,tolval});
    
    ang=4;
    % leave out 1:1206
    [tmp1(ang,:,:),tmp2] = STRFestimate({stim_c4(:,1207:end)},{resp_c4(N,1207:end)},100,{0,{30,0:29},'DirectFit'},{sparseval,tolval});
    
    ang=5;
    % leave out 1:1462
    [tmp1(ang,:,:),tmp2] = STRFestimate({stim_c5(:,1463:end)},{resp_c5(N,1463:end)},100,{0,{30,0:29},'DirectFit'},{sparseval,tolval});
    
    
    save(['strf_subj1_ch_',num2str(N)],'tmp1');
    figure;
    cc = max(abs(tmp1(:)));
    for j=1:5
        subplot(2,3,j)
        cc = max2(abs(tmp1(j,:,:)));
        imagesc(squeeze(tmp1(j,:,:)),[-cc cc]); colormap(jet); title(['Channel',num2str(N),'Angle',num2str(j)]);
        
    end
    k=[];
    k=strcat(num2str(N),'sub1.jpg');
    saveas(1,k); close all;
end

%%

for N=1:64 % loop channel
    k=strcat('strf_subj1_ch_',num2str(N));
    k=strcat(k,'.mat');
    tmp1=load([k]);
    tmp1=tmp1.tmp1;
    for i =1:5 % loop strf
        for j=1:5 % loop resp
            if j==1
                %[presp,crr]  = predictresp(squeeze(tmp1(i,:,:)),stim_c1,resp_c1);
                [presp]  = predictresp(squeeze(tmp1(i,:,:)),stim_c1);
                [r,p] = corr(presp',resp_c1(N,:)');
                save(['R_resp1_ch_',num2str(N),'_strf',num2str(i)],'r');
                save(['P_resp1_ch_',num2str(N),'_strf',num2str(i)],'p');
                % save(['crrnum_resp1_ch_',num2str(N),'_strf',num2str(i)],'crr');
                save(['Pred_resp1_ch_',num2str(N),'_strf',num2str(i)],'presp');
            end
            
            if j==2
                %[presp,crr]  = predictresp(squeeze(tmp1(i,:,:)),stim_c2,resp_c2);
                [presp]  = predictresp(squeeze(tmp1(i,:,:)),stim_c2);
                [r,p] = corr(presp',resp_c2(N,:)');
                save(['R_resp2_ch_',num2str(N),'_strf',num2str(i)],'r');
                save(['P_resp2_ch_',num2str(N),'_strf',num2str(i)],'p');
                % save(['crrnum_resp2_ch_',num2str(N),'_strf',num2str(i)],'crr');
                save(['Pred_resp2_ch_',num2str(N),'_strf',num2str(i)],'presp');
            end
            
            if j==3
                %[presp,crr]  = predictresp(squeeze(tmp1(i,:,:)),stim_c3,resp_c3);
                [presp]  = predictresp(squeeze(tmp1(i,:,:)),stim_c3);
                [r,p] = corr(presp',resp_c3(N,:)');
                save(['R_resp3_ch_',num2str(N),'_strf',num2str(i)],'r');
                save(['P_resp3_ch_',num2str(N),'_strf',num2str(i)],'p');
                % save(['crrnum_resp3_ch_',num2str(N),'_strf',num2str(i)],'crr');
                save(['Pred_resp3_ch_',num2str(N),'_strf',num2str(i)],'presp');
            end
            
            
            if j==4
                %[presp,crr]  = predictresp(squeeze(tmp1(i,:,:)),stim_c4,resp_c4);
                [presp]  = predictresp(squeeze(tmp1(i,:,:)),stim_c4);
                [r,p] = corr(presp',resp_c4(N,:)');
                save(['R_resp4_ch_',num2str(N),'_strf',num2str(i)],'r');
                save(['P_resp4_ch_',num2str(N),'_strf',num2str(i)],'p');
                %  save(['crrnum_resp4_ch_',num2str(N),'_strf',num2str(i)],'crr');
                save(['Pred_resp4_ch_',num2str(N),'_strf',num2str(i)],'presp');
            end
            
            if j==5
                % [presp,crr]  = predictresp(squeeze(tmp1(i,:,:)),stim_c5,resp_c5);
                [presp]  = predictresp(squeeze(tmp1(i,:,:)),stim_c5);
                [r,p] = corr(presp',resp_c5(N,:)');
                save(['R_resp5_ch_',num2str(N),'_strf',num2str(i)],'r');
                save(['P_resp5_ch_',num2str(N),'_strf',num2str(i)],'p');
                %  save(['crrnum_resp5_ch_',num2str(N),'_strf',num2str(i)],'crr');
                save(['Pred_resp5_ch_',num2str(N),'_strf',num2str(i)],'presp');
            end
        end
    end
end
%imst(tmp1);drawnow;
%end