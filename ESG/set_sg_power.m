function set_sg_power(sg,power)
% set_sg_power.m
%
% PURPOSE: To set power for the signal generator.
% PARAMETERS: sg - GPIB instrument object already initialized and open.
%             power - Output power.

% if (nargin~=2)
%     help set_sg_power
%     error('Wrong input arguments');
% end
% 
% sg=instrfind('Tag',['sg' num2str(sg_num)]);

fprintf(sg,[':POW ' num2str(power)]); % Setup output power.
%fprintf(sg,':OUTP ON'); % Turn RF on.
%fprintf(sg,':POW:ATT:AUTO ON') % % Turn Off attenuator hold.



