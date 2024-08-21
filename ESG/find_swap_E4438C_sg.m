function find_swap_E4438C_sg()

sg_i=0;
instr=instrfind;
for i=1:length(instr)
    instr_curr=instr(i).Tag;
    instr_type=instr_curr(1:2);
    if (strcmpi(instr_type,'sg'))
        sg_i=sg_i+1;
        sg(sg_i)=instr(i);
        fprintf(sg(sg_i),'*IDN?');
        gg=scanstr(sg(sg_i),',');
        g{sg_i}=gg{2};
    end
end


if (sg_i==0) 
    error('There are no signal generators');
end

if (sg_i>1) % If we have more than one generator
    for ii=1:(sg_i)
        if (strcmpi(g{ii},'E4438C') && ii~=1)
            sg(1).Tag=sg(ii).Tag;
            sg(ii).Tag='sg1';
        else 
            break;
        end
    end
else
    if(strcmpi(g,'E4438C'))
        error('E4438C not existent')    
    end    
end
