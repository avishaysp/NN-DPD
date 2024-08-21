function ESG = initiate_ESG_with_ipaddr(addr)
%% start interface with SG
%ESG_ip_addr = tcpip('132.68.61.205',5025);
ESG = tcpip(addr{1}, str2double(addr{2}));
ESG.OutputBufferSize = 600000;
ESG.ByteOrder = 'bigEndian';
ESG.Timeout = 10.0;
fopen(ESG);
query(ESG,'*IDN?')
fprintf(ESG,':OUTPut:STATe OFF')
fprintf(ESG,':SOURce:RADio:ARB:STATe OFF')
fprintf(ESG,':OUTPut:MODulation:STATe OFF')

end

