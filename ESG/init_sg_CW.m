function init_sg_MTONE(sg_num,freq)
% init_sg_CW.m
%
% PURPOSE: To initialize the signal generator to produce a CW signal.
% PARAMETERS: sg_num - number of sig. gen. .
%             freq - Center frequency.

if (nargin~=2)
    help init_sg_CW
    error('Wrong input arguments');
end

sg=instrfind('Tag',['sg' num2str(sg_num)]);
set(sg,'Timeout',300); % ???

fprintf(sg,'*IDN?'); % Get ID.
gg=scanstr(sg,',');
model=gg{2};
disp(model);

set_sg_off(sg_num); % Turn RF power off.
fprintf(sg,[':FREQ:FIX ' num2str(freq) ' MHz']); % Setup frequency.
fprintf(sg,':POW -50'); % Setup output power.

fprintf(sg,':OUTPut:MODulation OFF');
