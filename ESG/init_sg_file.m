function init_sg_file(sg_num,sig_name,sample_clock)
% init_sg_file.m
%
% PURPOSE: To initialize the signal generator to produce an arbitrary waveform.
%          It is assumed that the waveform file is already in the generator's non-volatile memory.
% PARAMETERS: sg_num - GPIB address for signal generator.
%                       sig_name - Name of arbitrary file to load;
%                       sample_clock - Sampling clock in MHz.

if (nargin~=3)
    help init_sg_file;
    error('Wrong input arguments');
end

sig_name

sg=instrfind('Tag',['sg' num2str(sg_num)]);
set(sg,'Timeout',300); % ???

%sig_name='802_11a_ADS_pre';

fprintf(sg,'*RST'); % Reset instrument.
fprintf(sg,'*IDN?'); % Get ID.
gg=scanstr(sg,',');
model=gg{2};
disp(model);

if(~strcmpi(model,'E4438C'))
    error('Can be used only with E4438C signal generator');
end

set_sg_off(sg_num); % Turn RF power off.
fprintf(sg,[':FREQ:FIX 2450 MHz']); % Setup frequency.
fprintf(sg,':POW -50'); % Setup output power.

fprintf(sg,':MMEM:CAT? ''WFM1'''); % Get file list for volatile memory.
fi=scanstr(sg);
fi_n=fi(3:3:end); % get only filenames.
ind_evdo=find(strcmpi(fi_n,['"' sig_name])==1);

if(isempty(ind_evdo))
    fprintf(sg,':MMEM:CAT? ''NVWFM'''); % Get file list for non-volatile memory.
    fi=scanstr(sg);
    fi_n=fi(3:3:end); % get only filenames.
    ind_evdo=find(strcmpi(fi_n,['"' sig_name])==1);
    if(isempty(ind_evdo))
        error('Signal is not available in generator volatile or non-volatile memory');
    end
    fprintf(sg,[':MEM:COPY "NVWFM:' sig_name '","WFM1:' sig_name '"']); % If it is in nv then copy it to v.
    fprintf(sg,'*OPC?') % Wait till the signal is constructed.
    read='';
    while(isempty(read))
        read=fscanf(sg,'%g');
    end
end

fprintf(sg,[':RADio:ARB:WAVeform "WFM1:' sig_name '"']); % Select the waveform.
fprintf(sg,[':RADio:ARB:SCLock:RATE ' num2str(sample_clock*1e6)]); % Set sample clock rate.
fprintf(sg,':RADio:ARB:HCRest ON'); % Turn on high crest mode.

fprintf(sg,':RADio:ARB ON'); % Turn on arbitraty generator.
fprintf(sg,':OUTPut:MODulation ON'); % Make sure modulation is on.

% fprintf(sg,'*OPC?') % Wait till the signal is constructed.
% read='';
% while(isempty(read))
%     read=fscanf(sg,'%g');
% end