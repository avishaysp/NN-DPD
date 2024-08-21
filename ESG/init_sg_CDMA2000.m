function init_sg_CDMA2000(sg_num,freq,num_carriers)
% init_sg.m
%
% PURPOSE: To initialize the signal generator to produce multicarrier CDMA2000-SR1.
% PARAMETERS: addr - GPIB adrress of instrument.
%             freq - Center frequency.
%             num_carriers - Number of CDMA2000-SR1 carriers to generate.

if (nargin~=3)
    help init_sg_CDMA2000
    error('Wrong input arguments');
end

sg=instrfind('Tag',['sg' num2str(sg_num)]);
set(sg,'Timeout',300); % ???
fprintf(sg,'*IDN?'); % Get ID.
%disp(fscanf(sg,'%*[^,], %[^,], %*[^,], %*[^,]'))
gg=scanstr(sg,',');
model=gg{2};
disp(model);

if ( strcmpi(model,'E4438C') || strcmpi(model,'ESG-D3000B') )   
    fprintf(sg,':DIAG:INFO:OPT?')
    str=scanstr(sg,',','%s');
    ok_flag=0;
    for i=1:length(str)
        if ( strcmpi(str{i},'401') || strcmpi(str{i},'101') )
                ok_flag=1;
                break;
        end
    end
    if (ok_flag==0)
        error('The signal generator doesn''t have CDMA2000 option');
    end
else
    error('Unknown signal generator')
end

set_sg_off(sg_num); % Turn RF power off.
fprintf(sg,[':FREQ:FIX ' num2str(freq) ' MHz']); % Setup frequency.
fprintf(sg,':POW -50'); % Setup output power.

% For IS-95
v=[1:num_carriers];
offset_car_vec=((v-(num_carriers+1)/2)*1.25e6); % Compute carrier offsets around center freq.

% Check current state of instrument to determine if signal has to be rebuilt.
fprintf(sg,':RADio:CDMA2000:ARB:LINK:FORWard:SETup?');
mode=fscanf(sg,'%s');
fprintf(sg,':RADio:CDMA2000:ARB:LINK:FORWard:SETup:MCARrier:TABLe:NCARriers?')
ncar=fscanf(sg,'%d');
fprintf(sg,':RADio:CDMA2000:ARB?');
stat=fscanf(sg,'%d');

if(strcmpi(mode,'MCAR') && ncar==num_carriers && stat) % We have multicarrier CDMA mode.
    for i=1:num_carriers % Check each carrier to see if it has the right offset
        fprintf(sg,[':RADio:CDMA2000:ARB:LINK:FORWard:SETup:MCARrier:TABLe? ' num2str(i)]); % Get the table row and see if it is what we need.
        row=scanstr(sg,',');
        if(row{2}~=offset_car_vec(i) || ~strcmpi(row{1},'S19CHAN'))
            break;
        end
    end
    if (i==num_carriers) % If we done all the carriers then the signal is ok.
        return;
    end
end

fprintf(sg,':RADio:CDMA2000:ARB:LINK:FORWard:SETup MCAR')
fprintf(sg,[':RADio:CDMA2000:ARB:LINK:FORWard:SETup:MCAR:TABL INIT,S19CHAN,' num2str(offset_car_vec(1)) ',0']);    
if(num_carriers>1)
    for i=2:num_carriers
        fprintf(sg,[':RADio:CDMA2000:ARB:LINK:FORWard:SETup:MCAR:TABL APP,S19CHAN,' num2str(offset_car_vec(i)) ',0']);
    end
end

fprintf(sg,':RADio:CDMA2000:ARB ON');

fprintf(sg,'*OPC?') % Wait till the signal is constructed.
read='';
while(isempty(read))
    read=fscanf(sg,'%g');
end
fprintf(sg,':OUTPut:MODulation ON');

