addpath('C:\Users\vlsi2\Documents\MATLAB\iqtools')
%% open communication with the MXG
ESG_addr = {'132.66.48.3','5025'};
ESG = initiate_ESG_with_ipaddr(ESG_addr);
%% Generate signal
sigBW = 4e6;
fs = 30e6;
modType = 'QAM16';
filterBeta = 0.25;
numSymbols = 4000;
iqdata = IQsigGen(sigBW, fs, modType, filterBeta, numSymbols);
iqdata = iqdata.';
numZeros = 3000 ;
iqdata = [zeros(1,numZeros) iqdata];
fRF = 500e6;
pRF = 9; %dBm

%%  write to signal generator
set_sg_freq(ESG, fRF) % center frequency
set_sg_power(ESG, pRF)   % output power

ESG_load_IQ(ESG, iqdata, fs*1e-6);
% set_sg_marker_Nimrod(ESG,1,'MATLAB_WFM.bin',1,round(length(iqdata)*0.75/100)*100,1)
set_sg_marker(ESG,1,'MATLAB_WFM.bin',1,round(length(iqdata)*0.75/100)*100,1)

%% open communication with the VSA
fclose(VSA);
VSA = open_vsa_visa();

%% Read IQ signal
trace = 3;
% fprintf(VSA, ':INPut:ANALOG:RANGe:AUTO');
Y = readVSA_visa_IQ(VSA, trace);
t = readVSA_time(VSA,trace);
EVM = db(str2double( ...
    query(VSA, ':TRACe4:DATA:TABLe:VALue? 1'))/100);

%% SYNC
close all;
outSig = Y/rms(Y);
fs_out = 1/(t(2)-t(1));

min_fs    = min(fs,fs_out);
iqdata_rs = resample(iqdata,double(round(min_fs)),double(round(fs)));
outSig_rs = resample(outSig,double(round(min_fs)),double(round(fs_out)));

[inAlign, outAlign, linear_gain_1] = Mini_Alignment(iqdata_rs.', outSig_rs.');

t_new = (0:1:length(inAlign)-1)/min_fs;

plot(t_new*1e6,real([inAlign outAlign])); grid on; ylabel 'Real';xlabel 't (\musec)'
% figure
% plot(t_new*1e6,imag([inAlign outAlign])); grid on; ylabel 'Imag';xlabel 't (\musec)'
%% PA modeling

% use stored wifi signals - not the calculated vectors above
%% Erez
linear_gain_1 = 1;
load('lotery_27_11.mat');
if isrow(sig_in)
    Input_signal = sig_in';
else
    Input_signal = sig_in;
end

if isrow(sig_out)
    Output_signal = sig_out';
    
else
    Output_signal = sig_out;
end

%% Past recorded signals of oursetup
load('captured_aligned_signals.mat');
Input_signal = inAlign;
Output_signal = outAlign;

%% Recorded signal
linear_gain_1=1;
Input_signal = inAlign;
Output_signal = outAlign * linear_gain_1;
% calc error vector
% Error magnitude
impedance = 50; % #FIXME
v_error = (Input_signal).*linear_gain_1 - (Output_signal); 
v_error_rms = 10*log10(rms(v_error).^2 / impedance * 1000);

%% PA Model - calc (G)MP coeffs
f = linspace(-64000000, 64000000, length(Input_signal)); %to check

K_PA = 3;    % PA fitting model - nonlinear 
M_PA = 3;     % PA fitting model - memory 
X_Matrix_PA = Build_Signal_Matrix(Input_signal, K_PA, M_PA);
PA_coef = Estimate_DPD_coeffs(X_Matrix_PA,Output_signal./linear_gain_1);

% Output_no_DPD_estimated
Output_no_DPD_estimated = X_Matrix_PA * PA_coef;
close all
if 1 % Plot AM/AM - AM/PM core
    figure; subplot(2,1,1); hold on; grid on;
    % 1. AM/AM RAW data 
    plot(10.*log10(linear_gain_1.*(abs(Input_signal)).^2/50*1000), ...
            20*log10( abs(Output_signal./linear_gain_1) ./ abs(Input_signal) ),'or','MarkerSize',2);
    % 2. AM/AM Fitted data
    plot(10.*log10(linear_gain_1 .* abs(Input_signal).^2/50*1000), ...
                20*log10(abs(Output_no_DPD_estimated)./abs(Input_signal)),'kx','MarkerSize',5);   
    axis([ -10 20 -10 5]);
    title('AM2AM Results'); xlabel('Pout [dBm]'); ylabel('AM/AM [dB]');

    subplot(2,1,2);  hold on; grid on;
    % 1. AM/PM RAW data 
    plot(10.*log10(linear_gain_1.*(abs(Input_signal)).^2/50*1000), ...
             ( angle(Output_signal./linear_gain_1 ./ Input_signal ) )./pi*180,'or','MarkerSize',2);        
    % 2. AM/PM Fitted data
    plot(10.*log10(linear_gain_1 .* abs(Input_signal).^2/50*1000), ...
            (angle(Output_no_DPD_estimated./Input_signal))./pi*180,'kx','MarkerSize',5);        
    xlabel('Pout [dBm]'); ylabel('AM/PM [Degrees]'); 
    title('AM2PM Results'); axis([-10 20 -50 50]);    
end
% model ranking
rank = normalizedSquaredDifferenceLoss(Output_no_DPD_estimated,Output_signal ./ linear_gain_1);
disp(['Ranking of the PA model: ', num2str(rank)]);
%% figures for the PA model frequency domain

figure;
plot(smooth(db(fftshift(fft(Input_signal))),500)); 
hold on;

plot(smooth(db(fftshift(fft(Output_signal))),500));
plot(smooth(db(fftshift(fft(Output_no_DPD_estimated))),500));
grid;
legend('Input Signal', 'Output Signal', 'Output without DPD Estimated');

NMSE_model = db(rms(Output_signal/rms(Output_signal)-Output_no_DPD_estimated/rms(Output_no_DPD_estimated)))

%% Calculate DPD coefficients
close all;
K_DPD = 7;
M_DPD = 7;
Y_Matrix = Build_Signal_Matrix(Output_signal./linear_gain_1, K_DPD, M_DPD);
DPD_coef_Inverse_Vector = Estimate_DPD_coeffs(Y_Matrix,Input_signal);

level = 0; % dB
X_Matrix = Build_Signal_Matrix(Input_signal .* (10.^(level./20)), K_DPD, M_DPD); 
Input_signal_DPD = X_Matrix*DPD_coef_Inverse_Vector;



% Model ranking
rank = normalizedSquaredDifferenceLoss(Input_signal_DPD, Input_signal);

% Display the rank
disp(['Model Rank: ', num2str(rank)]);



% plot(real([Input_signal, Output_signal ./ linear_gain_1, Input_signal_DPD])); grid on; ylabel 'Real';xlabel 't (\musec)'
% legend('Input Signal', 'Output Signal / Linear Gain 1', 'Input after DPD');

PAPR = db(max(abs(Input_signal_DPD))/rms(Input_signal_DPD))


%% figures for the DPD model frequency domain and AM/AM AM/PM

Y_after_DPD = Y_Matrix * DPD_coef_Inverse_Vector;
figure;
subplot(3,1,1);
plot(smooth(db(fftshift(fft(Output_signal))),500)); 
hold on;

plot(smooth(db(fftshift(fft(Input_signal))),500));
plot(smooth(db(fftshift(fft(Y_after_DPD))),500));
grid;
legend('Output Signal', 'Input Signal', 'Y after DPD');
 
subplot(3,1,2); hold on; grid on;
% 1. AM/AM RAW data 
plot(10.*log10(abs(Output_signal).^2/50*1000), 20*log10( abs(Input_signal) ./ abs(Output_signal) ),'or','MarkerSize',2);
% 2. AM/AM Fitted data
plot(10.*log10(abs(Output_signal).^2/50*1000), 20*log10(abs(Y_after_DPD ) ./ abs(Output_signal)),'kx','MarkerSize',5);   
axis([ -10 20 -10 5]);
title('AM2AM Results'); xlabel('Pout [dBm]'); ylabel('AM/AM [dB]');

subplot(3,1,3);  hold on; grid on;

% 1. AM/PM RAW data 
plot(10.*log10(abs(Output_signal).^2/50*1000), (angle(Input_signal ./ Output_signal))./pi*180,'or','MarkerSize',2);        
% 2. AM/PM Fitted data
plot(10.*log10(abs(Output_signal).^2/50*1000), (angle(Y_after_DPD ./ Output_signal))./pi*180,'kx','MarkerSize',5);
xlabel('Pout [dBm]'); ylabel('AM/PM [Degrees]'); 
title(['AM2PM Results']); axis([-10 20 -50 50]);    

%% pass X through DPD, then PA. second try :)

X_Matrix_before_dpd = Build_Signal_Matrix(Input_signal, K_DPD, M_DPD);
input_signal_after_dpd = X_Matrix_before_dpd * DPD_coef_Inverse_Vector;

input_dpded_Matrix = Build_Signal_Matrix(input_signal_after_dpd ,K_PA, M_PA);
input_dpded_then_paed = input_dpded_Matrix * PA_coef;

%% figures for the Full Flow - frequency domain and AM/AM AM/PM

figure;
plot(smooth(db(fftshift(fft(Input_signal))),500)); hold on;
plot(smooth(db(fftshift(fft(Output_signal))),500));
plot(smooth(db(fftshift(fft(input_dpded_then_paed))),500));
grid;
legend('Input Signal', 'Output Signal', 'input DPDed then PAed');
%%
subplot(3,1,2); hold on; grid on;
% 1. AM/AM RAW data 
plot(10.*log10(abs(Output_signal).^2/50*1000), 20*log10( abs(Input_signal) ./ abs(Output_signal) ),'or','MarkerSize',2);
% 2. AM/AM Fitted data
plot(10.*log10(abs(Output_signal).^2/50*1000), 20*log10(abs(Y_after_DPD ) ./ abs(Output_signal)),'kx','MarkerSize',5);   
axis([ -10 20 -10 5]);
title('AM2AM Results'); xlabel('Pout [dBm]'); ylabel('AM/AM [dB]');

subplot(3,1,3);  hold on; grid on;

% 1. AM/PM RAW data 
plot(10.*log10(abs(Output_signal).^2/50*1000), (angle(Input_signal ./ Output_signal))./pi*180,'or','MarkerSize',2);        
% 2. AM/PM Fitted data
plot(10.*log10(abs(Output_signal).^2/50*1000), (angle(Y_after_DPD ./ Output_signal))./pi*180,'kx','MarkerSize',5);
xlabel('Pout [dBm]'); ylabel('AM/PM [Degrees]'); 
title(['AM2PM Results']); axis([-10 20 -50 50]);    




%% take a signal x and pass it through dpd filter and then a pa model

X_Matrix_DPDed_TO_PA = Build_Signal_Matrix(Input_signal_DPD,model,'PA'); % here we took the dpded signal in order to pass it throgh the pa coeffiscient
Output_with_DPD_estimated = X_Matrix_DPDed_TO_PA *PA_coef;




if 1 % Plot AM/AM - AM/PM core
    figure; subplot(2,1,1); hold on; grid on;
    % 1. AM/AM RAW data 
    plot(10.*log10(linear_gain_1.*(abs(Input_signal)).^2/50*1000), ...
            20*log10( abs(Output_signal./linear_gain_1) ./ abs(Input_signal) ),'or','MarkerSize',2);
    % 2. AM/AM Fitted data
    plot(10.*log10(linear_gain_1 .* abs(Input_signal).^2/50*1000), ...
                20*log10(abs(Output_with_DPD_estimated )./abs(Input_signal)),'kx','MarkerSize',5);   
    axis([ -10 20 -10 5]);
    title('AM2AM Results'); xlabel('Pout [dBm]'); ylabel('AM/AM [dB]');

    subplot(2,1,2);  hold on; grid on;

    % 1. AM/PM RAW data 
    plot(10.*log10(linear_gain_1.*(abs(Input_signal)).^2/50*1000), ...
             ( angle(Output_signal./linear_gain_1 ./ Input_signal ) )./pi*180,'or','MarkerSize',2);        
    % 2. AM/PM Fitted data
    plot(10.*log10(linear_gain_1 .* abs(Input_signal).^2/50*1000), ...
            (angle(Output_with_DPD_estimated ./Input_signal))./pi*180,'kx','MarkerSize',5);        
    xlabel('Pout [dBm]'); ylabel('AM/PM [Degrees]'); 
    title(['AM2PM Results']); axis([-10 20 -50 50]);    
end



% model ranking
rank = normalizedSquaredDifferenceLoss(Output_with_DPD_estimated,Output_signal ./ linear_gain_1);
disp(['Ranking of the PA model: ', num2str(rank)]);
%% figures for the PA model frequency domain
figure;
plot(smooth(db(fftshift(fft(Input_signal))),100));
hold on;
plot(smooth(db(fftshift(fft(Output_signal))),100));
plot(smooth(db(fftshift(fft(Output_with_DPD_estimated))),100));
grid;

NMSE_model = db(rms(Output_signal/rms(Output_signal)-Output_with_DPD_estimated/rms(Output_with_DPD_estimated)))

%% Send DPDed sinal to generator
set_sg_power(ESG,pRF + level);
ESG_load_IQ(ESG, Input_signal_DPD.', fs*1e-6);
set_sg_marker(ESG,1,'MATLAB_WFM.bin',1,round(length(Input_signal_DPD)*0.75/100)*100,1);
%% write predistorted signal to generator
