function init_sg_MTONE(sg_num,freq,num_tones,spacing)
% init_sg_MTONE.m
%
% PURPOSE: To initialize the signal generator to produce a multitone signal. It is assumed to be peak phased.
% PARAMETERS: addr - GPIB adrress of instrument.
%             freq - Center frequency.(MHz)
%             num_tones - Number of multitones to generate.
%             spacing - Frequency spacing between tones.(MHz)
%example :init_sg_MTONE(1,freq,20,800e3); % Set signal generator to 20 CW tones at 800 KHz apart.

if (nargin~=4)
    help init_sg_MTONE
    error('Wrong input arguments');
end

sg=instrfind('Tag',['sg' num2str(sg_num)]);
set(sg,'Timeout',300); % ???

fprintf(sg,'*IDN?'); % Get ID.
gg=scanstr(sg,',');
model=gg{2};
disp(model);

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