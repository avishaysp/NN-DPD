function init_sg_CDMA(sg_num,freq,num_carriers)
% init_sg.m
%
% PURPOSE: To initialize the signal generator to produce multicarrier CDMA2000-SR1.
% PARAMETERS: addr - GPIB adrress of instrument.
%             freq - Center frequency.
%             num_carriers - Number of CDMA2000-SR1 carriers to generate.

if (nargin~=3)
    help init_sg_CDMA
    error('Wrong input arguments');
end

sg=instrfind('Tag',['sg' num2str(sg_num)]);
set(sg,'Timeout',300); % ???
fprintf(sg,'*IDN?'); % Get ID.
gg=scanstr(sg,',');
model=gg{2};
disp(model);

if ( strcmpi(model,'E4438C') || strcmpi(model,'ESG-D3000B') )   
    fprintf(sg,':DIAG:INFO:OPT?')
    str=scanstr(sg,',','%s');
    ok_flag=0;
    for i=1:length(str)
        if ( strcmpi(str{i},'401') || strcmpi(str{i},'UN5') )
                ok_flag=1;
                break;
        end
    end
    if (ok_flag==0)
        error('The signal generator doesn''t have CDMA option');
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
fprintf(sg,':RADio:CDMA:ARB:SETup?');
mode=fscanf(sg,'%s');
fprintf(sg,':RADio:CDMA:ARB:STATe?');
stat=fscanf(sg,'%d');

if(strcmpi(mode,'MCAR') && stat==1)
    fprintf(sg,':RADio:CDMA:ARB:SETup:MCAR:TABL?'); % Get the table and see if it is what we need.
    tabl=scanstr(sg,','); % Formated scan with delimeter ','.
    if(mod(length(tabl),num_carriers)==0)
        if(length(tabl)>0)
            if(tabl{1}=='FWD9') % If we have the right number of carriers and the right type.
                offs=[tabl{3:4:end}]; % Get the carrier offsets.
                if(length(offs)==num_carriers)
                    if(offs==offset_car_vec) % If offsets match then we are ok and we have nothing to do.
                        return;
                    end
                end
            end
        end
    end
end

st=[];
for i=v
    st=[st ',FWD9,"",' num2str(offset_car_vec(i)) ',0'];
end
st=st(2:end);
 
fprintf(sg,':RADio:CDMA:ARB:SETup MCAR')
fprintf(sg,[':RADio:CDMA:ARB:SETup:MCAR:TABL ' st])
fprintf(sg,':RADio:CDMA:ARB:STATe ON');

fprintf(sg,'*OPC?') % Wait till the signal is constructed.
read='';
while(isempty(read))
    read=fscanf(sg,'%g');
end
fprintf(sg,':OUTPut:MODulation ON');
