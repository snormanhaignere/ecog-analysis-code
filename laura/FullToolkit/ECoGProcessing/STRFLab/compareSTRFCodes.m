% need to load a STRF first
% load('/Users/LauraLong/Desktop/out_Local_HighGamma_025_LIJ109.mat')
% tval = calcSpeechTVal_out(out);

channels = [];
gap = [];
plotvals = 1;
useSTRFLab = 1;
numchan = 5;

[a,b] = sort(tval,'descend');
bestchan = b(1:numchan);

strfs1 = calcSTRFs_out(out,bestchan,gap,plotvals,0);
strfs2 = calcSTRFs_out(out,bestchan,gap,plotvals,1);

%% 

figure; plot(mean(strfs1.r_values),'o-'); hold on; plot(mean(strfs2.r_values),'*-');
legend('James','STRFLab');

for i = 1:length(bestchan)
    figure,
    subplot(121)
    imagesc(squeeze(mean(strfs1.strf(:,:,:,i)))); axis xy
    title(['James''s STRF: Corr ' num2str(mean(strfs1.r_values(:,i)))]);
    subplot(122)
    imagesc(squeeze(mean(strfs2.strf(:,:,:,i)))); axis xy
    title(['STRFLab: Corr ' num2str(mean(strfs2.r_values(:,i)))]);
    colormap jet;
    suptitle(num2str(bestchan(i)))
end