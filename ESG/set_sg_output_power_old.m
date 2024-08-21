function set_sg_output_power(sg,pm,opower,gain,offsets)
% set_sg_output_power.m
%
% PURPOSE: To find input power for the signal generators that will generate
%          the desired output power. Supports any number of signal generators.
% PARAMETERS: sg - Array of GPIB instrument objects already initialized and open for signal generators.
%             pm - GPIB instrument object already initialized and open for power meter.
%             opower - Required output power in dBm.
%             gain - System gain. (Has to be within 3-4 dB).
%             offsets - Vector of offsets between generators. All relative
%                       to 1st generator. Should be same size as number of
%                       generators. First offset should always be 0.

if (nargin~=5)
    help set_sg_output_power
    error('Wrong input arguments');
end

if(length(sg)~=length(offsets))
    error('Number of generators does not match number of provided offsets');
end

num_gen=length(sg);

fprintf(sg,':POW:ATT:AUTO ON'); % Turn on attenuator hold.

stat=get_sg_status(sg);
if(stat)
    p=get_pm(pm);
    if(abs(p-opower)<2) % If we are already near the required ouput power and the sg is on then just go forward.
        power=get_sg_power(sg(1));
    end
end
if(~exist('power','var'))
    power=opower-gain-2-10*log10(num_gen);
    for i=1:num_gen
        set_sg_power(sg(i),power-offsets(i)); % Set power on all generators.
        set_sg_on(sg(i)); % Turn them on.
    end
end

p=-100;
c=0;
while(1)
    while(p<opower)
        for i=1:num_gen
            set_sg_power(sg(i),power-offsets(i)); % Set power.
        end
        p=get_pm(pm); % Measure output power.
        if(p<(opower-2)) % Some conditions to make it go fast.
            inc=1;
        elseif(p<(opower-0.2))
            inc=0.1;
        else inc=0.02;
        end
        power=power+inc;
    end
    if((p-opower)>0.05)
        c=c+1; % Check to see if we are not stuck in a loop due to the signal generator's attenuator problems.
        if(c>1) 
            fprintf(sg,':POW:ATT:AUTO OFF'); % Turn off attenuator hold.
        end
        while(p>opower)
            for i=1:num_gen
                set_sg_power(sg(i),power-offsets(i)); % Set power.
            end
            p=get_pm(pm); % Measure output power.
            power=power-0.02;
        end
    else
        break;
    end
end      