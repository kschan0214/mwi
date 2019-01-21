%% [output_thresholded,U,S,V] = svd3D_thres(input,mask,thres_component)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 15 october 2018
% Date last modified:
%
%
function [output_thresholded,U,S,V] = svd3D_thres(input,mask,thres_component)

% get matrix size
matrixSize = size(input);

% compute 3D svd
[U,S,V,~] = svd3D(input,mask);

% thresholding
S(thres_component:end,thres_component:end) = 0;
% compute threshold martix
input_thresholded_2D = U*S*V';

% restore matrix to input size
% input_thresholded_2D = reshape(input_thresholded_2D,[size(input_thresholded_2D,1) matrixSize(4:end)] );
tmp= zeros(matrixSize(1:3));
output_thresholded=zeros([matrixSize(1:3) prod(matrixSize(4:end))]);
for k= 1:size(input_thresholded_2D,2)
    tmp(mask==1) = input_thresholded_2D(:,k);
    output_thresholded(:,:,:,k) = tmp;
end
output_thresholded = reshape(output_thresholded,matrixSize);

output_U=zeros([matrixSize(1:3) prod(matrixSize(4:end))]);
for k= 1:size(input_thresholded_2D,2)
    tmp(mask==1) = U(:,k);
    output_U(:,:,:,k) = tmp;
end


U = output_U;



% for kfa=1:nfa
%     for kte=1:nt
%         tmp1(mask_thres==1)=input_thresholded_2D(:,kte,(kfa));
%         output_thresholded(:,:,:,kte,(kfa))=tmp1;
%     end
% end


end