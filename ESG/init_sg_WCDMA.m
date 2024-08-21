function init_sg_WCDMA(sg_num,freq,num_carriers,offset,tm,clipping)
% init_sg_WCDMA.m
%
% PURPOSE: To initialize the signal generator to produce multicarrier WCDMA. 
% PARAMETERS: sg_num - 
%             freq - Center frequency.
%             num_carriers - Number of WCDMA carriers to generate.
%             offset - offset frequency between carriers [in Hz]. Default=5 MHz.
%             tm - Test model to use. If not specified defaults to TM1D64. Values are 1:TM1D64, 4:TM4.
%             clipping - Perform clipping for 8.5 dB PAR signal. For the moment defined only for TM1D64.
%                        Default=0, otherwise specify clipping value (50% for 8.5 dB PAR).

if (nargin<3 || nargin>6)
    help init_sg_WCDMA
    error('Wrong input arguments');
end

if (nargin<5)
    tm=1;
end
if (nargin<4)
    offset=5e6;
end

carr_type={'TM1D64','','','TM4'};

if (nargin<6)
    clipping=0;
end

sg=instrfind('Tag',['sg' num2str(sg_num)]);
set(sg,'Timeout',300); % ???
fprintf(sg,'*IDN?'); % Get ID.
gg=scanstr(sg,',');
model=gg{2};
disp(model);

if (nargin==6 && ~strcmpi(model,'E4438C') && num_carriers~=1)
    warning('This setup does not support clipping');
    clipping=0;
end

%%%%%%%    checking if the signal generator have the opsion  %%%%%
if ( strcmpi(model,'E4438C') || strcmpi(model,'ESG-D3000B') )   
    fprintf(sg,':DIAG:INFO:OPT?')
    str=scanstr(sg,',','%s');
    ok_flag=0;
    for i=1:length(str)
        if ( strcmpi(str{i},'400') || strcmpi(str{i},'100') )
            ok_flag=1;
            break;
        end
    end
    if (ok_flag==0)
        error('The signal generator doesn''t have WCDMA option');
    end
else
    error('Unknown signal generator')
end

set_sg_off(sg_num); % Turn RF power off.
fprintf(sg,[':FREQ:FIX ' num2str(freq) ' MHz']); % Setup frequency.
fprintf(sg,':POW -50'); % Setup output power.

if (num_carriers==1 && clipping)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(sg,':RADio:WCDMA:TGPP:ARB:LINK:DOWN:SETup?'); % Get current setup.
    mode=fscanf(sg,'%s');
    fprintf(sg,':RADio:WCDMa:TGPP:ARB:CLIPping?'); % Get clipping level for single carrier.
    clipp=fscanf(sg,'%g');
        
    if(strcmpi(mode,'TM1D64') && clipp==clipping) % IF the mode is single carrier with the right clipping then do nothing.
        return;
    end
    
    fprintf(sg,':RADio:WCDMA:TGPP:ARB:LINK:DOWN:SETup TM1D64');
    fprintf(sg,[':RADio:WCDMa:TGPP:ARB:CLIPping ' num2str(clipping)]);
    fprintf(sg,':RADio:WCDMa:TGPP:ARB ON');
    fprintf(sg,'*OPC?') % Wait till the signal is constructed.
    read='';
    while(isempty(read))
        read=fscanf(sg,'%g');
        pause(1);
    end
else
    v=[1:num_carriers];
    offset_car_vec=((v-(num_carriers+1)/2)*offset); % Compute carrier offsets around center freq.
    
    % Check current state of instrument to determine if signal has to be rebuilt.
    fprintf(sg,':RADio:WCDMA:TGPP:ARB:LINK:DOWN:SETup?');
    mode=fscanf(sg,'%s');
    fprintf(sg,':RADio:WCDMA:TGPP:ARB:LINK:DOWN:SETup:MCARrier:TABLe:NCARriers?')
    ncar=fscanf(sg,'%d');
    fprintf(sg,':RADio:WCDMA:TGPP:ARB?');
    stat=fscanf(sg,'%d');
    
    diff=0;
    if(strcmpi(mode,'MCAR') && ncar==num_carriers && stat) % We have multicarrier WCDMA mode.
        for i=1:num_carriers % Check each carrier to see if it has the right offset
            fprintf(sg,[':RADio:WCDMA:TGPP:ARB:LINK:DOWN:SETup:MCARrier:TABLe? ' num2str(i)]); % Get the table row and see if it is what we need.
            row=scanstr(sg,',');
            if(row{2}~=offset_car_vec(i) || ~strcmpi(row{1},carr_type{tm}))
                diff=1;
            end
        end
        if (~diff) % If we done all the carriers then the signal is ok.
            return;
        end
    end
    
    fprintf(sg,':RADio:WCDMA:TGPP:ARB:LINK:DOWN:SETup MCAR') % Multicarrier mode.
    if(clipping)
        fprintf(sg,[':RADio:WCDMA:TGPP:ARB:LINK:DOWN:SETup:MCAR:TABL INIT,' carr_type{tm} ',' num2str(offset_car_vec(1)) ',0,0,0,0,' num2str(clipping) ',100']);    
    else
        fprintf(sg,[':RADio:WCDMA:TGPP:ARB:LINK:DOWN:SETup:MCAR:TABL INIT,' carr_type{tm} ',' num2str(offset_car_vec(1)) ',0']);
    end
    if(num_carriers>1)
        for i=2:num_carriers
            if(clipping)
                fprintf(sg,[':RADio:WCDMA:TGPP:ARB:LINK:DOWN:SETup:MCAR:TABL APP,' carr_type{tm} ',' num2str(offset_car_vec(i)) ',0,0,0,0,' num2str(clipping) ',100']);
            else
                fprintf(sg,[':RADio:WCDMA:TGPP:ARB:LINK:DOWN:SETup:MCAR:TABL APP,' carr_type{tm} ',' num2str(offset_car_vec(i)) ',0']);    
            end
        end
    end
    
    if(strcmpi(model,'E4438C'))
        fprintf(sg,':RADio:WCDMa:TGPP:ARB:LINK:DOWN:SETup:MCARrier:SCODe:AINCrement');   %auto-increment scramble codes
        fprintf(sg,':RADio:WCDMa:TGPP:ARB:LINK:DOWN:SETup:MCARrier:TOFFset:AINCrement'); %auto-increment timing offsets
    end
    %fprintf(sg,':RADio:WCDMa:TGPP:ARB:LINK:DOWN:SETup:TABL:APPLy');
    fprintf(sg,':RADio:WCDMa:TGPP:ARB ON');
    
    fprintf(sg,'*OPC?') % Wait till the signal is constructed.
    read='';
    while(isempty(read))
        read=fscanf(sg,'%g');
        pause(1);
    end
    %fprintf(sg,':OUTPut:MODulation ON');
end