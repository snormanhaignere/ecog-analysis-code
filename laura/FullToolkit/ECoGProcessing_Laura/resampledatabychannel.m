function [newdata] = resampledatabychannel(data,upsamp,downsamp)
% [newdata] = resampledatabychannel(data,upsamp,downsamp)
% data must be chan x trial x sample

origsize = size(data);
newdata = zeros(origsize(1),origsize(2),origsize(3)/downsamp*upsamp);

for i = 1:origsize(1)
    disp(i);
    thiselec = data(i,:,:);
    newelec = resample(squeeze(thiselec)',1,2);
    newelec = permute(newelec,[3 2 1]);
    
    newdata(i,:,:) = newelec;
    
end

end

