%% [output,residual] = polyfit_image3D(input,model,mask)
%
% Input
% --------------
% input         : 3D image
% model         : 4D image, first 3 dimensions have the same size as the
%                 input, the last dimension contains polyfit term(s)
% mask          : 3D image with the same size as 
%
% Output
% --------------
% output        : 3D fitted image
% residual      : 3D difference between the fitted result and the input
%
% Description: Fitting the input 3D image to the provided model
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 15 october 2018
% Date last modified:
%
%
function [output,residual] = polyfit_image3D(input,model,mask)
% if mask is not provided 
if nargin < 3
    mask = ones(size(input));
end

% get matrix size
matrixSize = size(input);

% check size of input
if matrixSize(1) ~= size(model,1) || matrixSize(2) ~= size(model,2) || matrixSize(3) ~= size(model,3)
    error('Size of input image does not match the size of model');
end
if matrixSize(1) ~= size(mask,1) || matrixSize(2) ~= size(mask,2) || matrixSize(3) ~= size(mask,3)
    error('Size of input image does not match the size of mask');
end

% get mask index
ind = mask == 1;
% reshape model from 4D to 2D - 1st dim: image voxels; 2nd dim: SVD terms
model = reshape(model,prod(matrixSize),size(model,4));
% apply mask to model
model = model(ind,:);
% apply mask to input image
input_1D = input(ind);

% compute coefficients 
b = pinv(model)*input_1D;

% use coefficients to compute fitted image
output_1D = model * b ;
output = zeros(matrixSize);
output(ind) = output_1D;

% compute residuals
residual = input-output;

end
