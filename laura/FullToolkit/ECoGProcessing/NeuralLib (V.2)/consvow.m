%  s2=CUSpeech;
%  s2 = set(s2,'PreStimSilence',.5);
%  s2 = set(s2,'PostStimSilence',.5);
function [] = consvow(s2,subject,out,channelnames, threshval,condition )
if exist('threshval') && strcmp(threshval,'default')
    threshval=10;
end

if ~exist('out') || isempty(out)
    load(['out_' subject '_2_15hz.mat']);
end

if ~exist('condition') || isempty(condition)
    condition='';
end

if isempty(dir(['Phn_' subject condition '.mat']))
    warning(' Phoneme data was not found... making it from out structure 2 to 15 hz ...');
    
    
    phns = SpeechEvents(s2,'Phonemes','list');
    temp=struct2cell(phns(1:40));
    phnsub=squeeze(temp(1,1,:));
    
    fname='resp';
    param=[0.3 0.6];
    
    phndata={};
    
    
    
    
    for i=1:length(phnsub)
        phoneme=phnsub(i)
        if strcmp(phoneme,'sp')
            phndata{i}=0;
            continue;
        end
        [phndata{i},~] = NeuralPhonemeRep(s2, out, phoneme, fname, param,'Single');%, feats,loc,gender);
    end
    
    ind=find(cellfun(@(x) strcmp(x,'sp'),phnsub));
    phnsub(ind)=[];
    phndata(ind)=[];
    
    save(['Phn_',subject condition],'phndata','phnsub')
else
    load(['Phn_' subject condition '.mat']);
    % check all channels
end

plosive = {'B','D','G','K','P','T'};

fricative = {'DH','F','S','SH','TH','V','Z'};
affricate={'JH','CH'};
nasal = {'M','N'};
vowel = {'AA','AO','OW','AH','UH','UW','IY','IH','EY','EH','AE','AW','AY',...
            'OY'};
% vowel = {'AA1','AE1','AH0','AH1','AO1','AW1','AY1','AY2','EH1','EH2','ER0',...
%     'ER1','EY1','EY2','IH0','IH1','IH2','IY0','IY1','OW0','OW1','OY1','UH1','UW1'};
%Ovowel={'AE1','AA1'};

semivowel={'HH','R','W','Y'}; %Approximant, lateral , HH
consonant = [plosive, fricative, nasal,affricate];
classes={'vowel','consonant'};



D = []; L = [];
for cnt1 = 1:length(phndata)
    
    if isempty(phndata{cnt1})
        continue
    end
    
    phoneme=phnsub(cnt1);
    
    flag=0;
    for k=1:length(classes)
        if any(ismember(eval(classes{k}),phoneme))
            j=k;
            flag=1;
            break;
        end
    end
    
    if flag==0
        display(phoneme);
        % continue
    end
    
    D = cat(3,D,phndata{cnt1}); % concatenate along 3rd dimension, in order of phonemes
    L = cat(1,L,repmat(j,size(phndata{cnt1},3),1));
end

display(['number of phonemes is: ' num2str(length(L))]);
D2 = reshape(D,[size(D,1) size(D,2)*size(D,3)]);
mean_D = mean(D2,2);
%    min_D = min(D2,[],2);
std_D  = std(D2,[],2);
D2 = D - repmat(mean_D,[1 size(D,2) size(D,3)]);
D2 = D2./repmat(std_D,[1 size(D,2) size(D,3)]);
D=D2;

for loop=1:2
    
    if exist('threshval')
        x=0;
        ind=[];
        for i=1:length(D)
            tmp=D(1:40,:,i);
            if max(tmp(:))>threshval
                max(tmp(:))
                x=x+1;
                ind=[ind i];
            end
        end
        D(:,:,ind)=[];
        L(ind)=[];
        display([num2str(length(ind)) ' phonemes are deleted becauese of artifact']);
    end
    
    D2 = reshape(D,[size(D,1) size(D,2)*size(D,3)]);
    mean_D = mean(D2,2);
    %    min_D = min(D2,[],2);
    std_D  = std(D2,[],2);
    D2 = D - repmat(mean_D,[1 size(D,2) size(D,3)]);
    D2 = D2./repmat(std_D,[1 size(D,2) size(D,3)]);
    D=D2;
    
end
%
%figure();
totalf=[];
D_normalized =D;

display('Calculating the f-ratio for each electrode ...');
for k=1:length(channelnames)
    %  subplot(10,7,k);
    f=[];
    electrodes=k;%[8 9 10 11 13 14 15 16 17 18 19 20 22 23 24 25 26 27 28 29 31 32 33 34 35 36 37 38 40 41 42 43];%20:240
    for i = 1:size(D_normalized,2);
        DP = squeeze(D_normalized(:,i,:));
        DP=DP(electrodes,:);
        [~, f(i) ,~] = DPrime(DP,L);
        %subplot(9,7,el);
        
        %[allf, f(i) ,fs(i)] = DPrime(DP,L,phns);
    end
    totalf(k,:)=f;
    
    %     plot(f);
    %     title([num2str(k) ' ' channelnames{k}]);
end
save(['totalf_' subject condition],'totalf');

t_p=[38 45 57];

t1=[];
t2=[];
for i=1:length(t_p)
    [t1(i,:),t2(i,:)]=sort(totalf(:,t_p(i)),'descend');
end

% FigHandle = figure;
% set(FigHandle, 'Position', [100, 100, 800, 300]);
% for i=1:length(t_p)
%
%     subplot(1,length(t_p),i);
%     topoplot(totalf(:,t_p(i)),'gtec62.locs');
%
%     caxis([0, max(t1(:))])
%     colorbar
%     h=title([' Time ' num2str(t_p(i)*10-300) ' ms']);
%     set(h,'FontSize',18)
%
% end

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 800, 600]);
for i=1:length(t_p)
    
    subplot(length(t_p),1,i);
    p=plot(linspace(-100,600,71),totalf(t2(i,1),20:90),'LineWidth',4);
    line([0 0],[0 50],'color','r')
    h=title([' Channel  ' num2str(t2(i,1))]);
    set(h,'FontSize',18)
    xlabel('Time (ms)');
    ylabel('F-ratio');
end

end
