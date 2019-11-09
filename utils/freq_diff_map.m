%% function fdm = freq_diff_map(img,t)
%
% Input
% --------------
% img           : complex-valued 4D image 
% t             : echo time
%
% Output
% --------------
% fdm           : Frequency difference map
%
% Description: Tendler and Bowtell MRM (2018) 104
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 25 July 2019
% Date last modified:
%
%
function fdm = freq_diff_map(img,t)

img_prime = bsxfun(@rdivide,img,img(:,:,:,1));
fdm = zeros(size(img_prime));
for k=1:length(t)
    fdm(:,:,:,k) = angle(img_prime(:,:,:,k) ./ img_prime(:,:,:,2).^(k-1)) / (2*pi*(t(k)-t(2)));
end
fdm = fdm(:,:,:,3:end);

end