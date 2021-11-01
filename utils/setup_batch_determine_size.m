%% [numBatch, nElement] = setup_batch_determine_size(mask,numBatch)
%
% Input
% --------------
% mask          : signal mask
% numBatch      : number of batches the data to be broken down
%
% Output
% --------------
% numBatch      : number of batches the data to be broken down (final)
% nElement      : number of elements in each batch
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 10 August 2021
% Date modified:
%
%
function [numBatch, nElement] = setup_batch_determine_size(mask,numBatch)

% find masked voxels
mask    = reshape(mask,numel(mask),1);
ind     = find(mask>0);

% repackaging for progress reporting and temporary storage
numMaskedVoxel = length(ind);

% check how many voxels per current batch
nElement = floor(numMaskedVoxel/numBatch);

% Scenario 1
% If the no. of voxels less than 200 or no. batches then do it in one go
if numMaskedVoxel < numBatch || numMaskedVoxel < 200
    numBatch = 1;
    % update nElement
    nElement = floor(numMaskedVoxel/numBatch);
    
    return
end

% Scenario 2
% avoid too many voxels per batch which waits long to update progress
while nElement > 500 && numBatch <=5000
    numBatch = numBatch + 50;
    nElement = floor(numMaskedVoxel/numBatch);
end

% avoid too few voxels per batch which reduces parallel computing efficiency
while nElement < 20 && numBatch ~= 1
    numBatch = round(numBatch/2);
    nElement = floor(numMaskedVoxel/numBatch);
end

% get final no. of voxles per batch
nElement = floor(numMaskedVoxel/numBatch);

end