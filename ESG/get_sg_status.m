function stat=get_sg_status(sg_num)
% get_sg_status.m
%
% PURPOSE: To indicate of signal generator is on or off.
% PARAMETERS: sg - GPIB instrument object already initialized and open.

if (nargin~=1)
    help get_sg_status
    error('Wrong input arguments');
end

sg=instrfind('Tag',['sg' num2str(sg_num)]);

fprintf(sg,':OUTP?'); % Get output status.
stat=fscanf(sg,'%d');



