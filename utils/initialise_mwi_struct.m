%% function output = function_name(input)
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
% Date created: 
% Date last modified:
%
%
function [algoParam, imgParam] = initialise_mwi_struct(varargin)

if isempty(varargin)
    method = 't2s';
else
    method = varargin{1};
end

algoParam.maxIter       = 1000;
algoParam.isParallel    = false;
algoParam.isWeighted    = false; % cost function weighted by abs(echo)/norm(echoes)
algoParam.DEBUG         = false;
algoParam.numMagn       = 1;       % mixed fitting

imgParam.te         = [];
imgParam.img        = [];
imgParam.mask       = [];
imgParam.fieldmap   = [];
    
%% GRE-VFA-T2* model
if ContainName(method,'vfa')
    
    algoParam.npulse = 50;        % no. of pulses for EPG-X
    
    imgParam.tr     = [];
    imgParam.fa     = [];
    imgParam.b1map  = [];
end


end
    

function bool = ContainName(name,string)
    bool= ~isempty(strfind(lower(name),string));
end