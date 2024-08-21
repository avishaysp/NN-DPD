function set_sg_on(ESG)
% set_sg.m
%
% PURPOSE: To turn RF power on for the signal generator.
% PARAMETERS: sg - GPIB instrument object already initialized and open.

% if (nargin~=1)
%     help set_sg_on
%     error('Wrong input arguments');
% end
% 
% sg=instrfind('Tag',['sg' num2str(sg_num)]);

fprintf(ESG,':OUTP ON'); % Turn RF on.
ESG.UserData='on';

