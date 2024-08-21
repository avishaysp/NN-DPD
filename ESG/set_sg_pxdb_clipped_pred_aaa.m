function set_sg_pxdb_clipped_pred(pxdb,o_name,sig_name,off_i_d,off_o_d,xnn_flag)
%  set_sg_pxdb_clipped_pred.m
%
% PURPOSE: To display the waveform at the desire PxdB.
% PARAMETERS: pxdb - the desire PxdB 
%                          o_name - name of Amp. for load the profile.
%                         sig_name - signal name
%                         off_i_d - offset table of the input detector
%                         off_o_d -offset table of the output detector
%                         xnn_flag - =1 with xnn ;  0 -without xnn

if (nargin~=6)
    help set_sg_pxdb_clipped_pred
    error('Wrong input arguments');
end


if xnn_flag==1
    load([o_name '\' o_name '_p1db_xnn.mat'])
else
    load([o_name '\' o_name '_p1db.mat'])
end

f_name=[o_name '\' sig_name];% '_pre_'  num2str(pxdb)];
sig_name_pre=[sig_name];%'_pre_'  num2str(pxdb)];
in=findstr(f_name,'.')
f_name(in)='_';
f_name
sig_name_pre
load(f_name)

clear ii qq
po_xdb
%po_xdb=27;

in=findstr(sig_name_pre,'.')
sig_name_pre(in)='_';
init_sg_file(1,sig_name_pre,80)                                  %
%ggg=robustfit(power_i(1:55),power_o(1:55));                  % Find linear fit for very small signal gain.
%pout_pxdb=interp1(power_o(end/2:end)-power_i(end/2:end),power_o(end/2:end),ggg(1)-pxdb);
set_sg_peak_output_power(1,po_xdb,off_i_d,off_o_d,30);
p_out=get_pm('pmo')
init_sg_file(1,sig_name,80)
set_sg_output_power(1,p_out-4,30);
[v1,v2]=get_sc_pkpk('CH1','CH2');
set_sc_scaling(v1,v2);
set_sg_output_power(1,p_out,30);
set_sc_trigger(1);
set_sc_ref_waveform(1,1)
init_sg_file(1,sig_name_pre,80)
set_sg_output_power(1,p_out,30);
set_sc_trigger(2);
%set_sc_ch_off(1)
