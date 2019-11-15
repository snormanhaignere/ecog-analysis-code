function sflag = writehtk(name,MAT,sampPeriod,paramKind) 
% writehtkf(name,MAT,sampPeriod,paramKind); 
%
% Simple function to write an HTK file. 
% HTK deafult (BIG endian) on the output. 
% MAT is the data, dim*samples
% samPeriod is the sampling frequency in Hz 

if ~exist('paramKind','var') || isempty(paramKind), 
    paramKind = 8971;
end
[P,nSamples]=size(MAT);
sampSize = P*4;
sampPeriod = sampPeriod*1E4;
%ff = fopen (name,'w','n');
ff = fopen(name,'w','b');  %default
if ff==-1, 
    sflag = 0; 
    return;
else
    sflag=1;
end
% write header
fwrite (ff,nSamples,'int');
fwrite (ff,sampPeriod,'int');
fwrite (ff,sampSize,'short');
fwrite (ff,paramKind,'short');
% determine amount of data to read and read 
fwrite (ff, MAT , 'float');
fclose (ff);
