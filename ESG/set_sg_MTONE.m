function set_sg_MTONE(sg_num,freq,num_tones,spacing)
% set_sg_MTONE.m
%
% PURPOSE: To initialize the signal generator to produce a multitone signal. It is assumed to be peak phased.
% PARAMETERS: sg - GPIB object already opened and initialized for signal generator.
%             freq - Center frequency [MHz].
%             num_tones - Number of multitones to generate.
%             spacing - Frequency spacing between tones.

if (nargin~=4)
    help set_sg_MTONE
    error('Wrong input arguments');
end

sg=instrfind('Tag',['sg' num2str(sg_num)]);

set_sg_off(sg_num); % Turn RF power off.
fprintf(sg,'*RST');
fprintf(sg,[':FREQ:FIX ' num2str(freq) ' MHz']); % Setup frequency.
fprintf(sg,':POW -50'); % Setup output power.

fprintf(sg,[':RADio:MTONe:ARB:SETup:TABLe ' num2str(spacing) ',' num2str(num_tones)]);
fprintf(sg,':RADio:MTONe:ARB:SETup:TABLe:PHASe:INITialize FIXed');
fprintf(sg,':RADio:MTONe:ARB ON')

fprintf(sg,'*OPC?') % Wait till the signal is constructed.
read='';
while(isempty(read))
    read=fscanf(sg,'%g');
end
fprintf(sg,':OUTPut:MODulation ON');
