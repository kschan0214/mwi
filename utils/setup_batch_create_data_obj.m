%% data_obj = setup_batch_create_data_obj(data, mask, numBatch, nElement, fieldName, data_obj)
%
% Input
% --------------
% data          : data to be store in batches, at least 2D, can have singleton dimension
% mask          : signal mask
% numBatch      : number of batches the data to be broken down (final)
% nElement      : number of elements in each batch
% fieldName     : field name to be store in data_obj, string
% data_obj      : (optional) previous created data_obj
% 
%
% Output
% --------------
% data_obj      : data_obj with data stored
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 10 August 2021
% Date modified:
%
%
function data_obj = setup_batch_create_data_obj(data, mask, numBatch, nElement, fieldName, data_obj)

% initialise data obj
if nargin < 6 || isempty(data_obj)
    data_obj = [];
end

% get matrix size of the data
dim         = size(data);
if numel(dim) == 1
    dim = [dim, 1, 1];
elseif numel(dim) == 2
    dim = [dim, 1];
end
% get the number of voxels
numVoxel    = prod(dim(1:3));

% find masked voxels
mask	= reshape(mask,numVoxel,1);
ind     = find(mask>0);

% Step 1: convert the matrix of the data such that 1st dim: voxel; 2nd/3rd dim: same as 4th/5th original
% Step 2: masked out the data
if numel(dim) == 3      % 3D image
    data = reshape(data,numVoxel,1);    
    data = data(ind);
    
elseif numel(dim) == 4  % 4D image
    data = reshape(data,numVoxel,dim(4));
    data = data(ind,:);
    
elseif numel(dim) == 5  % 5D image
    data = reshape(data,numVoxel,dim(4),dim(5));
    data = data(ind,:,:);
    
end

% store data into batches
for kbat = 1:numBatch
    startInd = (kbat-1)*nElement+1;
    
    if kbat < numBatch
        endInd      = kbat*nElement;
    else
        endInd      = length(ind);
    end
    
    if numel(dim) == 3      % 3D image
        data_obj(kbat).(fieldName) = data(startInd:endInd);
        
    elseif numel(dim) == 4  % 4D image
        data_obj(kbat).(fieldName) = data(startInd:endInd,:);
        
    elseif numel(dim) == 5  % 5D image
        data_obj(kbat).(fieldName) = data(startInd:endInd,:,:);
        
    end
    
end


end