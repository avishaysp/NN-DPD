function set_sg_output_power(num_gen,opower,assumed_gain)
% set_sg_output_power.m
%
% PURPOSE: To find input power for the signal generators that will generate
%          the desired output power. Supports any number of signal generators.
% PARAMETERS: num_gen - the number of signal generators that are needed
%             opower - Required output power in dBm.

if (nargin~=2 && nargin~=3)
    help set_sg_output_power
    error('Wrong input arguments');
end

if(nargin==2)
    assumed_gain=80;
end

opower=opower-10*log10(num_gen)+0.02;
opower_orig=opower;

sg_i=1;
for g=1:30
    sgt=instrfind('Tag',['sg' num2str(sg_i)]);
    if( ~isempty(sgt) )
        sg(sg_i)=sgt;
        sg_i=sg_i+1;
    else
        break;
    end
end
if ( (sg_i-1) < num_gen)
    error('There are not enough signal generators');
end

pmi=instrfind('Tag','pmi');
pmo=instrfind('Tag','pmo');

%fprintf(pmi,'SENS:AVER:COUNT 32') % Set averaging filter for input power meter to 32 to increase measurement speed.

for i=1:num_gen
    fprintf(sg(i),':POW:ATT:AUTO ON'); % Turn off attenuator hold.
end

for i=1:num_gen
    stat=get_sg_status(i); % See if the generator is on.
    if(stat)
        power(i)=get_sg_power(i); % If so, then use the current power.
        set_sg_off(i); % Turn it off.
    else
        power(i)=opower-assumed_gain-10*log10(num_gen); % Assume maximum gain of 80 dB.
        set_sg_power(i,power(i)); % Set poweron generator.        
    end
end


lll=0;
while(1)
    err_flag=0;
    for g=1:num_gen % Loop over all generators
        if(err_flag) % If error flag is set then we have overpower condition and we have to leave this loop.
            break
        end
        set_sg_on(g); % Turn generator on.   
        p_old1=0; % Used to detect if we are stuck due to attenuator jump.
        p_old2=0;
        pi_prev=get_pm(pmi.tag);
        po_prev=get_pm(pmo.tag);
        pg_prev=power(g);
        for i=1:20 % Maximum Newton-Raphson iterations.
            set_sg_power(g,power(g)); % Set power.
            pause(0.5);
            p=get_pm(pmo.tag);
            d=p-opower;
            if(abs(d)<0.03) 
                break;
            end
            %%%% Sanity check - Test whether the gain change at input and output is consistent with the gain change in the generator.
            pi=get_pm(pmi.tag); % Get input power for error checking.
            gi=pi-pi_prev; % Gain delta for input.
            go=p-po_prev; % Gain delta for output.
            gg=power(g)-pg_prev; % Gain delta for generator.
            if( (abs(go-gg)>3 || abs(gi-gg)>5) && pi>-5 ) % Check gain change.
                set_sg_off(g); % Turn off generator.
                warning('Something is wrong - The gain of either input or output is not consistent with the generator');
                fprintf('gg:%+5.2f (%+5.2f : %+5.2f), gi: %+5.2f (%+5.2f : %+5.2f), go:%+5.2f (%+5.2f : %+5.2f)\n',gg,pg_prev,power(g),gi,pi_prev,pi,go,po_prev,p);
                while(1)
                    key=input('What to do? [Continue/Stop]','s');
                    switch key
                        case 's'
                            set_sg_off(g);
                            error('Program Stopped');
                        case 'c'
                            power(g)=power(g)-5; % Reduce power by 5 dB.
                            set_sg_power(g,power(g)); % Set the new power.
                            set_sg_on(g); % Turn the generator on.
                            err_flag=1;
                            break;
                    end
                end
            end
            pi_prev=pi;
            po_prev=p;
            pg_prev=power(g);
            
            power(g)=power(g)-d; % Newton-Raphson update.

            if( p > (opower_orig+2) && d<0) % Error trapping - Check if any generator has overpower of more than 2 dB
                set_sg_off(g); % Turn off the offending generator
                disp(['Something is wrong - Power level at output is too high: ' num2str(p) ' dBm.']); 
                fprintf('Powers: ');
                for kk=1:num_gen  % For debugging
                    fprintf('%5.2f dBm, ',power(kk));
                end
                fprintf('\n');
                disp('Waiting one minute and retrying');
                pause(60) % Wait one minute to recuperate
                err_flag=1; % Set error flag so that we may break from all inner loops
                power(g)=power(g)-5; % Reduce power in the offending generator by 5 dB so that when we restart it will not be too high.
                break; % Break from inner Newton-Raphson loop
            end
            if(abs(p-p_old2)<0.1) % We are probably stuck due to attenuator jump.
                fprintf(sg(g),':POW:ATT:AUTO OFF'); % Turn on attenuator hold.
            end
            p_old2=p_old1;
            p_old1=p;
        end
        if (i==20)
            warning('Maximum number of Newton-Raphson iterations reached without solution');
        end
        set_sg_off(g);
    end
    
    for i=1:num_gen % Turn on all generators.
        set_sg_on(i);
    end
    
    if(num_gen==1) 
        break;
    end
    pause(0.2);
    p=get_pm(pmo.tag);
    p=p-10*log10(num_gen);
    if(abs(p-opower_orig)>0.03)
        if(lll==0)
            lll=1;
            for i=1:num_gen
                set_sg_power(i,get_sg_power(i)-(p-opower))
            end
            pause(0.2)
            p=get_pm(pmo.tag);
            p=p-10*log10(num_gen);
            if(abs(p-opower_orig)>0.03)
                opower=opower_orig-(p-opower);
            else
                break;
            end
        else
            opower=opower_orig-(p-opower);
        end
    else
        break;
    end
    for i=1:num_gen % Turn off all generators.
        set_sg_off(i);
    end
end

fprintf(pmi,'SENS:AVER:COUNT:AUTO ON') % Set averaging filter for input power meter to original value.
fprintf(sg,':POW:ATT:AUTO ON') % % Turn Off attenuator hold.