function data_out = rereference(data_in,refold,refnew)
% REREFERENCE.M: re-reference neural data data
%
% INPUTS
%
%   data_in:        original data (electrodes x time)
%
%   refold:         channel index of old reference electrode(s)
%                   (for g.tec system, should be in the range of [1,..,64])
%                   can be scalar or vector
%                   optional; default is 64 (right earlobe for g.tec)
%
%   refnew:         channel index of target reference electrode(s)
%                   optional; default is [63 64] (linked earlobes for g.tec)
%
% OUTPUTS
%
%   data_out:       re-reference data
%
% by tasha
% in progress (HAS NOT BEEN TESTED)

if ~exist('refold', 'var') 
    refold = 64; 
end
if ~exist('refnew', 'var') 
    refnew = [63 64]; 
end

% remove old reference
channel = mean(data_in(refold,:),1);
Rmat = repmat(channel,[size(data_in,1) 1]); 
Rmat(refold,:) = 0; 
data_out = data_in + Rmat; 

% add new reference
channel = mean(data_out(refnew,:),1); 
Rmat = repmat(channel,[size(data_out,1) 1]); 
Rmat(refnew,:) = 0; 
data_out = data_out + Rmat; 

end