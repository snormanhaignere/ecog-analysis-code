load('/Users/LauraLong/Documents/Lab/ECoG Data/029_LIJ112/processed/B26/out_029_LIJ112_B26_htkraw.mat')
%%

close all; 

n = 5;

numchan = 1;

figure; 
hold on; 
subplot(numchan+3,1,1); 
plot(out(n).sound); 
title('sound'); 
subplot(numchan+3,1,2);
plot(out(n).audrecord);
title('sound on gtec');
subplot(numchan+3,1,3);
plot(out(n).reference);
title('reference');
for i = 1:numchan
    c = 64-i;
    subplot(numchan+3,1,i+3);
    plot(out(n).resp(c,:));
    title(['elec ' num2str(c)]);
end
suptitle(['emospeech trial ' num2str(n)]); 
hold off;

keyboard;
soundsc(zscore(out(n).reference),out(n).dataf);
keyboard;
soundsc(out(n).sound,out(n).soundf)
keyboard;
soundsc(zscore(out(n).reference),out(n).dataf);