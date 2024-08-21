function init_sg_15FA_CFR(sg_num)
% init_sg_15FA_CFR.m
%
% PURPOSE: To initialize the signal generator to produce a 15 carrier CFR source saved in nvram. 
%          It is assumed that the waveform file is already in the generator's non-volatile memory.
%          It is named: "chipx24_15FA_1frame". 
% PARAMETERS: sg_num - GPIB address for signal generator.

if (nargin~=1)
    help init_sg_15FA_CFR(addr);
    error('Wrong input arguments');
end

%find_swap_E4438C_sg; % Set E4438C to be sg1.

sg=instrfind('Tag',['sg' num2str(sg_num)]);
set(sg,'Timeout',300); % ???

fprintf(sg,'*RST'); % Reset instrument.
fprintf(sg,'*IDN?'); % Get ID.
gg=scanstr(sg,',');
model=gg{2};
disp(model);

if(~strcmpi(model,'E4438C'))
    error('Can be used only with E4438C signal generator');
end

set_sg_off(sg_num); % Turn RF power off.
fprintf(sg,[':FREQ:FIX 2140 MHz']); % Setup frequency.
fprintf(sg,':POW -50'); % Setup output power.

fprintf(sg,':MMEM:CAT? ''WFM1'''); % Get file list for volatile memory.
fi=scanstr(sg);
fi_n=fi(3:3:end); % get only filenames.
ind_evdo=find(strcmpi(fi_n,'"chipx24_15FA_1frame')==1);

if(isempty(ind_evdo))
    fprintf(sg,':MMEM:CAT? ''NVWFM'''); % Get file list for non-volatile memory.
    fi=scanstr(sg);
    fi_n=fi(3:3:end); % get only filenames.
    ind_evdo=find(strcmpi(fi_n,'"chipx24_15FA_1frame')==1);
    if(isempty(ind_evdo))
        error('Signal is not available in generator volatile or non-volatile memory');
    end
    fprintf(sg,':MEM:COPY "NVWFM:chipx24_15FA_1frame","WFM1:chipx24_15FA_1frame"'); % If it is in nv then copy it to v.
    fprintf(sg,'*OPC?') % Wait till the signal is constructed.
    read='';
    while(isempty(read))
        read=fscanf(sg,'%g');
    end
end

fprintf(sg,':RADio:ARB:WAVeform "WFM1:chipx24_15FA_1frame"'); % Select the waveform.
fprintf(sg,':RADio:ARB:SCLock:RATE 29.4912e6'); % Set sample clock rate.
fprintf(sg,':RADio:ARB:HCRest ON'); % Turn on high crest mode.

fprintf(sg,':RADio:ARB ON'); % Turn on arbitraty generator.
fprintf(sg,':OUTPut:MODulation ON'); % Make sure modulation is on.

% fprintf(sg,'*OPC?') % Wait till the signal is constructed.
% read='';
% while(isempty(read))
%     read=fscanf(sg,'%g');
% end