%% theta = AngleBetweenV1MapAndB0(v1,b0dir)
%
% Input
% --------------
% v1            : 4D fibre orientation map in vector form
% b0dir         : 1D vector of B0 direction
%
% Output
% --------------
% theta         : 3D angle map
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 20 March 2019
% Date last modified:
%
%
function theta = AngleBetweenV1MapAndB0(v1,b0dir)

% replicate B0 direction to all voxels
b0dirmap = permute(repmat(b0dir(:),1,size(v1,1),size(v1,2),size(v1,3)),[2 3 4 1]);
% compute angle between B0 direction and fibre orientation
theta = atan2(vecnorm(cross(v1,b0dirmap),2,4), dot(v1,b0dirmap,4));

end