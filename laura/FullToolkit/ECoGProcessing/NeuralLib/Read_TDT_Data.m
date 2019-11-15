DataPath = 'C:/Users/Laura/Documents/Columbia/Lab/ECoG Data/CUECoGRoDu/TDT_Data'
ReadFile = '/B3'

[data] = SEV2mat([DataPath ReadFile],'VERBOSE',0);
aud = double(data.Aud_.data(1,:));
fs_aud = double(data.Aud_.fs);
y = double(data.xWav.data);
fs_ECoG = double(data.xWav.fs);