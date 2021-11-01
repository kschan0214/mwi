%% [scaleFactor, img] = mwi_image_normalisation(img, scaleMask)
%
% Input
% --------------
% img           : image to be normalised
% scaleMask     : signal mask
%
% Output
% --------------
% scaleFactor   : scaling factor
% img           : normalised images
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 16 Nov 2020
% Date modified:
%
%
function [scaleFactor, img] = mwi_image_normalisation(img, scaleMask)

% reshape to 4D image
img_max = reshape(abs(img),[size(img,1),size(img,2),size(img,3),size(img,4)*size(img,5)]);
img_max = max(img_max,[],4);

% compute a signal mask if it is not provided
if nargin < 2
    scaleMask = img_max/prctile(img_max(:),95) > 0.15;
end

% determine the scaling factor
scaleFactor = norm(img_max(scaleMask>0)) / sqrt(length(find(scaleMask>0)));

img = img / scaleFactor;

end