function init_sg_1xEVDO(sg_num)
% init_sg_1xEV-DO.m
%
% PURPOSE: To initialize the signal generator to produce a 1xEV-DO 75% signal. 
%          It is assumed that the waveform file is already in the generator's non-volatile memory.
%          It is named: "EV-DO_75-WFM0". 
% PARAMETERS: addr - GPIB address for signal generator.

if (nargin~=1)
    help init_sg_1xEVDO
    error('Wrong input arguments');
end

sg=instrfind('Tag',['sg' num2str(sg_num)]);
set(sg,'Timeout',300); % ???

fprintf(sg,'*RST'); % Reset instrument.
fprintf(sg,'*IDN?'); % Get ID.
gg=scanstr(sg,',');
model=gg{2};
disp(model);

if ( strcmpi(model,'E4438C'))   
    fprintf(sg,':DIAG:INFO:OPT?')
    str=scanstr(sg,',','%s');
    ok_flag=0;
    for i=1:length(str)
        if ( strcmpi(str{i},'404'))
                ok_flag=1;
                break;
        end
    end
    if (ok_flag==0)
        error('The signal generator doesn''t have 1xEVDO option');
    end
else
    error('The signal generator is not E4438C')
end

set_sg_off(sg_num); % Turn RF power off.
fprintf(sg,[':FREQ:FIX 2140 MHz']); % Setup frequency.
fprintf(sg,':POW -50'); % Setup output power.

fprintf(sg,':MMEM:CAT? ''WFM1'''); % Get file list for volatile memory.
fi=scanstr(sg);
fi_n=fi(3:3:end); % get only filenames.
ind_evdo=find(strcmpi(fi_n,'"EV-DO_75-WFM0')==1);

if(isempty(ind_evdo))
    fprintf(sg,':MMEM:CAT? ''NVWFM'''); % Get file list for non-volatile memory.
    fi=scanstr(sg);
    fi_n=fi(3:3:end); % get only filenames.
    ind_evdo=find(strcmpi(fi_n,'"EV-DO_75-WFM0')==1);
    if(isempty(ind_evdo))
        error('Signal is not available in generator volatile or non-volatile memory');
    end
    fprintf(sg,':MEM:COPY "NVWFM:EV-DO_75-WFM0","WFM1:EV-DO_75-WFM0"'); % If it is in nv then copy it to v.
    fprintf(sg,'*OPC?') % Wait till the signal is constructed.
    read='';
    while(isempty(read))
        read=fscanf(sg,'%g');
    end
end

fprintf(sg,':RADio:ARB:WAVeform "WFM1:EV-DO_75-WFM0"'); % Select the waveform.
fprintf(sg,':RADio:ARB:SCLock:RATE 4.9152e6'); % Set sample clock rate.
fprintf(sg,':RADio:ARB:HCRest ON'); % Turn on high crest mode.

fprintf(sg,':RADio:ARB ON'); % Turn on arbitraty generator.
fprintf(sg,':OUTPut:MODulation ON'); % Make sure modulation is on.

% fprintf(sg,'*OPC?') % Wait till the signal is constructed.
% read='';
% while(isempty(read))
%     read=fscanf(sg,'%g');
% end
