function set_sg_marker(sg_num,marker,filename,start_m,end_m,periodic)
% set_sg_marker1.m
%
% PURPOSE: To set marker1 on the arbtraty eaveform to a periodic marker.
%          Marker will be off between start_m and end_m and will be on
%          between end_m and 2*start_m and so on until the end of the
%          waveform.
% PARAMETERS: sg_num - Number of signal generator.
%             marker - Marker Number 
%             filename - Name of arbitrary waveform signal that you wish to apply the marker to.
%             start_m - The start sample number in the waveform to apply the marker.
%             end_m - The end sample number in the waveform to apply the marker.
%             periodic - Define the marker to be periodic. (1=periodic). Default.

if (nargin<5)
    help set_sg_marker
    error('Wrong input arguments');
end

if(nargin==5)
    periodic=1;
end

% sg=instrfind('Tag',['sg' num2str(sg_num)]);
sg = sg_num;

fprintf(sg,':MMEMory:CATalog? "WFM1"'); % Get file list from volatile memory
fff=scanstr(sg,',');
i_fff=find(strcmpi(fff,['"' filename]));
fff2=fff{i_fff+2};
fff2=fff2(1:end-1);
n=str2num(fff2)/4;

fprintf(sg,['RADio:ARB:MARKer:CLEar:ALL "' filename '",' num2str(marker)])
skip_m=start_m;
on_m=end_m-start_m;
nn=n/(on_m+skip_m);
if(periodic)
    for ii=1:floor(nn)
        fprintf(sg,['RAD:ARB:MARK:SET "' filename '",' num2str(marker) ',' num2str(start_m) ',' num2str(end_m) ',0']);
        pause(0.1) 
        start_m=end_m+skip_m;
        end_m=start_m+on_m;
    end
    fprintf(sg,['RAD:ARB:MPOL:MARK' num2str(marker) ' NEG']);
else
    fprintf(sg,['RAD:ARB:MARK:SET "' filename '",' num2str(marker) ',' num2str(start_m) ',' num2str(end_m) ',0']);
end
