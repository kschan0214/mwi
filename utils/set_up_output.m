%% [output_dir,output_filename, temp_fn, identifier] = set_up_output(imgPara,output_filename,temp_prefix)
%
% Input
% --------------
% imgPara       : structure array contains all image data
% Output setting
% ==============
%   .output_dir      : directory to store final results (default:
%                      '/current_directory/mwi_results/')
%   .output_filename : output filename in text string (default:
%                      'mwi_results.mat')
%   .identifier      : temporary file identifier, a 8-digit code (optional)
%
% Output
% --------------
% output_filename : full output filename with extension
% temp_filename   : full temporary filename
% identifier      : 8-digit code (in string)
%
% Description: setup output related operations
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 16 Nov 2020
% Date modified:
%
%
function [output_dir,output_filename, temp_filename, identifier] = set_up_output(imgPara,output_filename,temp_prefix)
rng('shuffle');

% check if user specified output directory
if isfield(imgPara,'output_dir')
    output_dir = imgPara.output_dir;
else
    output_dir    = fullfile(pwd,'mwi_results');
end
if ~exist(output_dir,'dir')
    mkdir(output_dir)
end

% check if user specified output filename
if isfield(imgPara,'output_filename')
    output_filename = imgPara.output_filename;
    [~,output_filename,~] = fileparts(output_filename);
    output_filename = [output_filename '.mat'];
end
output_filename = fullfile(output_dir,output_filename);

% check if user specified idendifier for previous progress
if isfield(imgPara,'identifier')
    % get temporary file identifier if provided
    identifier  = imgPara.identifier;
    temp_dir    =  fullfile(output_dir,['temporary_' identifier]);
    temp_prefix = fullfile(temp_dir,temp_prefix);

    % check if the temporary file exists
    if ~exist([temp_prefix identifier '.mat'],'file')
        error('Cannot detect the temporary file with the provided identifier. Please enter a valid identifier or remove the input identifier');
    end
else
    % create a new identifier if not provided
    identifier = [];
    % make sure the identifier is unique
    while isempty(identifier) || exist([temp_prefix identifier '.mat'],'file')
        % identifier = num2str(randi(9));
        identifier = datestr(datetime('now'),'yymmddHHMMSSFFF');
        % for k = 2:8; identifier = [identifier,num2str(randi(9))]; end
        % identifier = num2str(identifier);
    end
    temp_dir    =  fullfile(output_dir,['temporary_' identifier]);
end
temp_filename = [temp_prefix identifier '.mat'];
% temp_filename = [temp_prefix identifier];

end