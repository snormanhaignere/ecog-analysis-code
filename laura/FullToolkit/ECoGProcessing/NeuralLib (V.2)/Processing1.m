%% Bahar 2015
% Generating phonetic transcription and Fratio for all channels

phns = SpeechEvents(s2,'Phonemes','list');


temp=struct2cell(phns(1:45));
phnsub=squeeze(temp(1,1,:));
phnsub(find(strcmp(phnsub,'sp')))=[];

fname='resp';
param=[0.3 0.6];

phndata={};


plosive = {'B','D','G','K','P','T'};

fricative = {'DH','F','S','SH','TH','V','Z'};
affricate={'JH','CH'};
nasal = {'M','N'};
vowel = {'AA1','AE1','AH0','AH1','AO1','AW1','AY1','AY2','EH1','EH2','ER0',...
    'ER1','EY1','EY2','IH0','IH1','IH2','IY0','IY1','OW0','OW1','OY1','UH1','UW1'};
Ovowel={'AE1','AA1'};

semivowel={'HH','R','W','Y'}; %Approximant, lateral , HH
consonant = [plosive, fricative, nasal,affricate];
classes={'vowel','consonant'};

phnsub=[plosive,fricative,affricate,semivowel,nasal,vowel];

for i=1:length(phnsub)
    phoneme=phnsub(i)
    [phndata{i},~] = NeuralPhonemeRep(s2, out, phoneme, fname, param,'Single');%, feats,loc,gender);
end

ind=[];
for i=1:length(phnsub)
    if length(phnsub{i})==1
        continue
    end
    
    c=~cellfun('isempty',strfind(phnsub,phnsub{i}(1:2)));
    c=find(c)
    if i==c(1)
        for j=2:length(c)
            ind=[ind c(j)];
            phndata{i}=cat(3,phndata{i},phndata{c(j)});
        end
    end
    phnsub{i}=phnsub{i}(1:2);
end
phndata(ind)=[];
phnsub(ind)=[];
%%
D = []; L = [];
for cnt1 = 1:length(phndata)
    
    if isempty(phndata{cnt1})
        continue
    end
    
    phoneme=phnsub(cnt1);
    
    %     flag=0;
    %     for k=1:length(classes)
    %         if any(ismember(eval(classes{k}),phoneme))
    %             j=k;
    %             flag=1;
    %             break;
    %         end
    %     end
    %
    %     if flag==0
    %         display(phoneme);
    %         % continue
    %     end
    
    D = cat(3,D,phndata{cnt1}); % concatenate along 3rd dimension, in order of phonemes
    L = cat(1,L,repmat(cnt1,size(phndata{cnt1},3),1));
end

%%
atlist = attribute2phoneme([],'list');
PL = zeros(length(unique(L)),size(D,3));
AL = zeros(length(atlist),size(D,3));
for cnt1 = 1:size(PL,2)
    PL(L(cnt1),cnt1)=1;
end
% find the attribs of this phoneme:
for cnt1 = 1:length(phnsub)
    atr = phoneme2attribute(phnsub{cnt1});
    atvec = zeros(1,length(atlist));
    for cnt2 = 1:length(atlist)
        if ~isempty(find(strcmpi(atr,atlist{cnt2}))), atvec(:,cnt2)=1;end
    end
    AL(:,L==cnt1)=repmat(atvec',[1 length(find(L==cnt1))]);
end
%%
f=[];
for j=1:62
electrodes=j;%[8 9 10 11 13 14 15 16 17 18 19 20 22 23 24 25 26 27 28 29 31 32 33 34 35 36 37 38 40 41 42 43];%20:240
Ltmp=AL(4,:);
for i = 1:90
    DP = squeeze(D(:,i,:));
    DP=DP(electrodes,:);
    [~, f(i) ,~] = DPrime(DP,Ltmp(:));
    %subplot(9,7,el);
    
    %[allf, f(i) ,fs(i)] = DPrime(DP,L,phns);
end
totalf(j,:)=f;
end
figure
plot(totalf')
line([30 30],[0 max(f)],'color','r','LineStyle','--');

[x t]=max(totalf(:,40));
figure
plot(totalf(t,:));
save('totalf.mat','totalf');
title(['Channel' num2str(t)]);
line([30 30],[0 max(totalf(:))],'color','r','LineStyle','--');