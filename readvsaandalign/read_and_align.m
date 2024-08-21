IP_addr = {'132.68.138.209','5026'};%IP_addr of VSA
trace = 5; %trace of VSA to read;

[XDelta, Y] = readVSA(IP_addr);
out_rms = rms(double(Y));
Output_signal = double(Y)/rms(double(Y));
Output_sample_rate = 1/Output_signal.XDelta*1e-6;

mutual_sample_rate   = min(Input_signal_SR,Output_sample_rate);
Input_signal         = resample(Input_signal,double(round(mutual_sample_rate*1e6)),double(round(signal.signal_sample_rate*1e6)));
Output_signal        = resample(Output_signal,double(round(mutual_sample_rate*1e6)),double(round(Output_sample_rate*1e6)));                  


[Input_signal,~,Output_signal] = Mini_Alignment(Input_signal, Output_signal);