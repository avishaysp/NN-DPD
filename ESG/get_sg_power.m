function power=get_sg_power(sg)
% get_sg_power.m
%
% PURPOSE: To read power setting on the signal generator.
% PARAMETERS: sg - GPIB instrument object already initialized and open.


% if (nargin~=1)
%     help get_sg_power
%     error('Wrong input arguments');
% end
% 
% sg=instrfind('Tag',['sg' num2str(sg_num)]);

fprintf(sg,':POW?'); % Get output power.
power=fscanf(sg,'%f');


