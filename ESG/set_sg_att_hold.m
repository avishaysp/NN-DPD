function set_sg_att_hold(sg_num,stat)
% set_sg.m
%
% PURPOSE: To set attenuator hold for the signal generator.
% PARAMETERS: sg - GPIB instrument object already initialized and open.
%             stat - 0:Att hold off. 1: Att hold on.

if (nargin~=2)
    help set_sg_att_hold
    error('Wrong input arguments');
end

sg=instrfind('Tag',['sg' num2str(sg_num)]);

if(stat)
    fprintf(sg,':POW:ATT:AUTO OFF'); % Turn on attenuator hold.
else
    fprintf(sg,':POW:ATT:AUTO ON'); % Turn on attenuator hold.
end



