%% function filename = check_unique_filename_seqeunce(filename)
%
% Input
% --------------
% filename      : filename to be checked if it is unique
%
% Output
% --------------
% filename      : a unique filename
%
% Description: to check if the input filename is unique or not, if not then
%              add a sequence number as the suffix of the name
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 16 Nov 2020
% Date modified:
%
%
function filename = check_unique_filename_seqeunce(filename)

[filepath,name,extension] = fileparts(filename);

if exist(filename,'file') == 2
    counter = 1;
    while exist(filename,'file') == 2
        suffix = ['_' num2str(counter)];
        filename = fullfile(filepath, [name suffix extension]);
        counter = counter + 1;
    end
end

end