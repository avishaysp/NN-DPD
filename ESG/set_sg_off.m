function set_sg_off(ESG)
% set_sg.m
%
% PURPOSE: To turn RF power off for the signal generator.
% PARAMETERS: sg - GPIB instrument object already initialized and open.

% if (nargin~=1)
%     help set_sg_off
%     error('Wrong input arguments');
% end
% 
% sg=instrfind('Tag',['sg' num2str(sg)]);

fprintf(ESG,':OUTP OFF'); % Turn RF off.
ESG.UserData='off';

