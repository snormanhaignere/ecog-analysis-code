function [compressout] = compressresp_out(out,whichcompressions,tanhfactor)
% [compressout] = compressresp_out(out,whichcompressions,tanhfactor)
%
% Compresses resp of an out structure
% Applies mapstd and tanh to all resp entries in out
% Takes tanhfactor or defaults to 10
% Takes whichcompressions which is 2x1 logical of [applyMapstd applyTanh]
%
% Written by Laura, NAPLab March 2017


fprintf('...compressing out structure');

if ~exist('tanhfactor','var') || isempty(tanhfactor)
    tanhfactor = 10;
    fprintf('\ntanhfactor defaulting to 10');
end
if ~exist('whichcompresssions','var') || isempty(whichcompressions)
    whichcompressions = [1 1];
    fprintf('\ndefaulting to mapstd and tanh compressions');
end



%% Compress each entry in out structure
compressout = out; % replace normout
for i = 1:length(out)
    thisresp = out(i).resp;
    if whichcompressions(1)
        thisresp = mapstd(thisresp);
    end
    if whichcompressions(2)
        thisresp = tanhfactor * tanh(thisresp / tanhfactor);
    end
    compressout(i).resp = thisresp;
end

fprintf('...done.\n');

end
