function ESGdn(sig_name,sample_clk)
% ESGdn.m
% 
% PURPOSE: To download a signal into the arbitrary  signal generator. 
% PARAMETERS: sig_name - Name of signal to load.
%                        sample_clk - The sample clock for the signal in MHz.
% for example   : ESGdn('.\WL#07\802_11a_ADS_pre_0_3',80)
root_dir=pwd;
path(root_dir,path);
tools_dir=strcat(root_dir,'\PSG-ESG Download Assistant');
path(tools_dir,path);

if(nargin~=2)
    help ESGdn
    error('Wrong input arguments');
end

%io=agt_newconnection('TCPIP','165.213.125.42');
io = agt_newconnection('gpib',0,16);
[status,status_descript,query_results]=agt_query(io,'*idn?')

load([sig_name '.mat'])
IQdata=ii'+j*qq';
clear ii qq;
disp('Done loading file');

slash=findstr(sig_name,'\');
if(~isempty(slash))
    if slash>=2
        l_slash=length(slash);
        begin=slash(l_slash)+1;
        sig_name=sig_name(begin:1:end);
    else    
        sig_name=sig_name(slash(end-1)+1:end);
    end
end
sig_name

%IQdata=IQdata/max(abs(IQdata))*0.7;
[status,status_descript]=agt_waveformload(io,IQdata,sig_name,sample_clk*1e6,'no_play')
[status, status_description] = agt_sendcommand(io,'SOURce:FREQuency 2450000000')
[status, status_description] = agt_sendcommand(io, 'POWer -30');
[status, status_description] = agt_sendcommand(io, [':RADio:ARB:WAVeform "WFM1:' sig_name '"']);
[status, status_description] = agt_sendcommand(io,[':RADio:ARB:SCLock:RATE ' num2str(sample_clk*1e6)]);
[status, status_description] = agt_sendcommand(io,':source:rad:arb:state on');
[status, status_description] = agt_sendcommand(io,':OUTPut:MODulation ON');

