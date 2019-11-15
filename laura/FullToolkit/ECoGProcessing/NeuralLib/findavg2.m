function [D_av,L_av,phnsub] = findavg2(s2,subject)

load(['out_' subject '_2_15hz']);

phns = get(s2,'Phonemes'); 
phnsub = cell(0,1); 
for cnt1 = 1:length(phns)
    phnstmp = struct2cell(phns{cnt1}); 
    phnstmp = squeeze(phnstmp(1,:,:)); 
    phnsub(end+1:end+length(phnstmp)) = phnstmp; 
end
phnsub = unique(phnsub); 

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

D = []; L = [];
for cnt1 = 1:length(phnsub)
    if isempty(phndata{cnt1}) || length(phndata{cnt1})==1
        continue
    end
    D = cat(3,D,phndata{cnt1}); % concatenate along 3rd dimension, in order of phonemes
    L = cat(1,L,repmat(cnt1,size(phndata{cnt1},3),1));
end

D2 = reshape(D,[size(D,1) size(D,2)*size(D,3)]);
mean_D = mean(D2,2);
min_D = min(D2,[],2);
std_D  = std(D2,[],2);
D2 = D - repmat(mean_D,[1 size(D,2) size(D,3)]);
D2 = D2./repmat(std_D,[1 size(D,2) size(D,3)]);
D=D2;

x=0;
ind=[];
for i=1:length(D)
    tmp=D(:,:,i);
    if max(tmp(:))>10
        max(tmp(:))
        x=x+1
        ind=[ind i];
    end
end
D(:,:,ind)=[];
L(ind)=[];

D_av=D;
L_av=L;



end

