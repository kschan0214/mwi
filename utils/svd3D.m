%% [U,S,V] = svd3D(input,mask)
%
% Input
% --------------
% input         : N-D image (N>3). Data after 4th dimension will be
%                 concatenated 
% mask          : 3D mask
%
% Output
% --------------
% U             : same as U in svd
% S             : same as S in svd
% V             : same as V in svd
% U_3D          : reshaped U, SVD components on the last dimension
%
% Description: perform 3D-based SVD
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 15 October 2018
% Date last modified:
%
%
function [U,S,V,U_3D] = svd3D(input,mask)

% get matrix size
matrixSize = size(input);

if length(matrixSize) < 3
    matrixSize(4) = 1;
end

% reshape data into 2D matrix - row: 3D image; col: other
input_2D = reshape(input,prod(matrixSize(1:3)),prod(matrixSize(4:end)));

% mask out non-target voxels
input_2D_masked = input_2D(mask(:)==1,:);

% econ svd
[U,S,V] = svd(input_2D_masked,0);

U_3D = zeros([matrixSize(1:3) prod(matrixSize(4:end))]);
tmp=zeros(matrixSize(1:3));

for ksvd=1:prod(matrixSize(4:end))
    tmp(mask==1)=U(:,ksvd);
    U_3D(:,:,:,ksvd)=tmp;
end

end