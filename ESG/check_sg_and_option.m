function [sg_e4438c,found]=check_sg_and_option(sg,signal)
% check_sg_and_option.m
%
% PURPOSE: To check if the signal generator has an option.
% PARAMETERS: sp - GPIB Instrument object already initialized and open.
%             signal - The signal we want to check for. If signal is not
%                      supplied, the function returns all the options for sg.
% RETURN: sg_e4438c - Boolean. 1 if the generator is E4438C.
%         found - Boolean. 1 if the option exists.
% for E4438C     : '400'-Wcdma    '401'-CDMA2000 and IS-95 CDMA    '404'-1xEV-DO
% for ESG-D3000B : '100'-Wcdma    '101'-CDMA2000                   'UN5'-IS-95

if (nargin~=2 && nargin~=1)
    help check_sg_and_option
    error('Wrong input arguments');
end


if nargin==1
    fprintf(sg,':DIAG:INFO:OPT?');
    option_list=scanstr(sg,',','%s')
else
    switch signal
        case {'wcdma','wcdma_8.5'}
            op={'400','100'};
        case 'is-95' 
            op={'401','UN5'};
        case '1xEV-DO'
            op={'404','404'}
        otherwise
            disp('Unknown option')
    end
    
    fprintf(sg,'*IDN?');
    gg=scanstr(sg,',');
    model=gg{2};
    
    if (strcmpi(model,'E4438C'))
        op_index=1;
        sg_e4438c=1;
    elseif (strcmpi(model,'ESG-D3000B'))
        op_index=2;
        sg_e4438c=0;
    else
        disp('Uknown generator');
        sg_e4438c=0;
        found=0;
        return;
    end
    
    fprintf(sg,':DIAG:INFO:OPT?');
    option_list=scanstr(sg,',','%s');
    found=0;
    for k=1:length(option_list)
        if ( strcmpi(option_list{k},op(op_index)) )
            found=1;
            break;
        end
    end
end