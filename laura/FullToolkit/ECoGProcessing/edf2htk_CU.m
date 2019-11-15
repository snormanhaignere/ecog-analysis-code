%% edf to htk files for CU recordings
% read .edf files
clear; close; clc;
%[hdr, rcd] = edfread('Passive20140819.edf'); 
[hdr, rcd] = edfread('KG.edf'); 
%%
%record = rcd(:,36000:210000);
%record = rcd(:,1200000:1370000);
%record = rcd(:,940000:1073000);
record = rcd(:,1073000:1200000);
% Audio channels 1
writehtk('./B04/htkFiles/Sync.htk', record(1,:),500);

% Brain recording channels
for cnt = 1:56
    filename = ['./B04/htkFiles/','Ch',num2str(cnt),'.htk'];
    writehtk(filename, record(cnt+1,:), 500);
end
filename = ['./B04/htkFiles/','ChAll.htk'];
writehtk(filename, record(2:57,:), 500);