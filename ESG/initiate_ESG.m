function ESG = initiate_ESG(addr)
%% start interface with SG
ESG = visa('agilent',addr);
ESG.OutputBufferSize = 100000;
ESG.ByteOrder = 'bigEndian';
ESG.Timeout = 30.0;
fopen(ESG);
idn = query(ESG, '*IDN?');
fprintf (idn);
end

